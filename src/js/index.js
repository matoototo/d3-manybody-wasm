import createModule from '../../dist/d3-manybody-wasm.js';

let moduleInstance = null;
let initializationPromise = null;

const STRIDE = 4;
const BYTES_PER_VALUE = Float32Array.BYTES_PER_ELEMENT;
const nodeBufferRegistry = new WeakMap();

function ensureModuleReady() {
    if (!moduleInstance) {
        throw new Error('WASM module not initialized. Please wait for initialization to complete.');
    }
}

function refreshEntryView(entry) {
    if (!entry || !entry.ptr || entry.count === 0) {
        entry.view = null;
        return null;
    }

    const desiredLength = entry.count * STRIDE;
    const heapBuffer = moduleInstance.HEAPF32.buffer;

    if (!entry.view || entry.view.buffer !== heapBuffer || entry.view.length !== desiredLength) {
        entry.view = new Float32Array(heapBuffer, entry.ptr, desiredLength);
    }

    return entry.view;
}

function retainNodeBuffer(nodes, count) {
    ensureModuleReady();

    let entry = nodeBufferRegistry.get(nodes);
    if (!entry) {
        entry = { ptr: 0, capacity: 0, view: null, refCount: 0, count: 0, lastUploadAlpha: Number.NaN, lastDownloadAlpha: Number.NaN };
        nodeBufferRegistry.set(nodes, entry);
    }

    const required = count * STRIDE;
    if (required > entry.capacity || entry.ptr === 0) {
        if (entry.ptr) {
            moduleInstance._free(entry.ptr);
        }
        if (required > 0) {
            entry.ptr = moduleInstance._malloc(required * BYTES_PER_VALUE);
            entry.capacity = required;
        } else {
            entry.ptr = 0;
            entry.capacity = 0;
        }
    }

    entry.count = count;
    entry.refCount += 1;
    refreshEntryView(entry);
    // Reset per-tick markers when (re)allocating/retaining
    entry.lastUploadAlpha = Number.NaN;
    entry.lastDownloadAlpha = Number.NaN;
    return entry;
}

function releaseNodeBuffer(nodes) {
    const entry = nodeBufferRegistry.get(nodes);
    if (!entry) return;
    entry.refCount -= 1;
    if (entry.refCount <= 0) {
        if (entry.ptr) {
            moduleInstance._free(entry.ptr);
        }
        nodeBufferRegistry.delete(nodes);
    }
}

function writeNodesToBuffer(nodes, view, count) {
    if (!view) return;
    for (let i = 0; i < count; ++i) {
        const node = nodes[i] || {};
        const base = i * STRIDE;
        view[base] = node.x ?? 0;
        view[base + 1] = node.y ?? 0;
        view[base + 2] = node.vx ?? 0;
        view[base + 3] = node.vy ?? 0;
    }
}

function writeBufferToNodes(view, nodes, count) {
    if (!view) return;
    for (let i = 0; i < count; ++i) {
        const node = nodes[i];
        if (!node) continue;
        const base = i * STRIDE;
        node.vx = view[base + 2];
        node.vy = view[base + 3];
    }
}

function initializeWasm() {
    if (!initializationPromise) {
        initializationPromise = createModule().then(module => {
            moduleInstance = module;
        });
    }
    return initializationPromise;
}

// Initialize WASM module immediately
initializeWasm();

function createAxisForce(createForceFunc, coordinateName) {
    return function axisForceFactory(coordinate) {
        ensureModuleReady();

        const force = createForceFunc(coordinate);
        let nodesRef = null;
        let nodeCount = 0;
        let bufferEntry = null;
        let bufferView = null;

        const syncNodesToBuffer = (alpha) => {
            if (!nodesRef || nodeCount === 0 || !bufferEntry) return;
            if (bufferEntry.lastUploadAlpha === alpha) return;
            bufferEntry.count = nodeCount;
            bufferView = refreshEntryView(bufferEntry);
            if (!bufferView) return;
            writeNodesToBuffer(nodesRef, bufferView, nodeCount);
            bufferEntry.lastUploadAlpha = alpha;
        };

        const syncBufferToNodes = (_alpha) => {
            if (!nodesRef || nodeCount === 0 || !bufferEntry) return;
            bufferView = refreshEntryView(bufferEntry);
            if (!bufferView) return;
            writeBufferToNodes(bufferView, nodesRef, nodeCount);
        };

        function releaseBuffer() {
            if (bufferEntry && nodesRef) {
                force.setNodeBuffer(0, 0);
                releaseNodeBuffer(nodesRef);
            }
            bufferEntry = null;
            bufferView = null;
        }

        function forceWrapper(alpha) {
            if (!nodesRef || nodeCount === 0 || !bufferEntry) return;
            syncNodesToBuffer(alpha);
            force.force(alpha);
            syncBufferToNodes(alpha);
        }

        forceWrapper.initialize = function (nodes) {
            releaseBuffer();
            nodesRef = nodes || null;
            nodeCount = nodesRef ? nodesRef.length : 0;
            bufferEntry = null;
            bufferView = null;

            force.setNodes(nodes);

            if (nodesRef && nodeCount > 0) {
                bufferEntry = retainNodeBuffer(nodesRef, nodeCount);
                if (bufferEntry.ptr) {
                    force.setNodeBuffer(bufferEntry.ptr, nodeCount);
                    bufferEntry.count = nodeCount;
                    bufferView = refreshEntryView(bufferEntry);
                }
            }

            if (!bufferEntry) {
                force.setNodeBuffer(0, 0);
            }

            return forceWrapper;
        };

        forceWrapper.dispose = function () {
            releaseBuffer();
            nodesRef = null;
            nodeCount = 0;
        };

        forceWrapper.strength = function (_) {
            if (arguments.length) {
                force.setStrength(_);
                return forceWrapper;
            }
            return force.getStrength();
        };

        forceWrapper[coordinateName] = function (_) {
            if (arguments.length) {
                force.setCoordinate(_);
                return forceWrapper;
            }
            return force.getCoordinate();
        };

        return forceWrapper;
    };
}

function createForceManyBody() {
    ensureModuleReady();

    const force = moduleInstance.createForceManyBody();
    let nodesRef = null;
    let nodeCount = 0;
    let bufferEntry = null;
    let bufferView = null;

    const syncNodesToBuffer = (alpha) => {
        if (!nodesRef || nodeCount === 0 || !bufferEntry) return;
        if (bufferEntry.lastUploadAlpha === alpha) return;
        bufferEntry.count = nodeCount;
        bufferView = refreshEntryView(bufferEntry);
        if (!bufferView) return;
        writeNodesToBuffer(nodesRef, bufferView, nodeCount);
        bufferEntry.lastUploadAlpha = alpha;
    };

    const syncBufferToNodes = (_alpha) => {
        if (!nodesRef || nodeCount === 0 || !bufferEntry) return;
        bufferView = refreshEntryView(bufferEntry);
        if (!bufferView) return;
        writeBufferToNodes(bufferView, nodesRef, nodeCount);
    };

    function releaseBuffer() {
        if (bufferEntry && nodesRef) {
            force.setNodeBuffer(0, 0);
            releaseNodeBuffer(nodesRef);
        }
        bufferEntry = null;
        bufferView = null;
    }

    function forceWrapper(alpha) {
        if (!nodesRef || nodeCount === 0 || !bufferEntry) return;
        syncNodesToBuffer(alpha);
        force.force(alpha);
        syncBufferToNodes(alpha);
    }

    forceWrapper.initialize = function (nodes) {
        releaseBuffer();
        nodesRef = nodes || null;
        nodeCount = nodesRef ? nodesRef.length : 0;
        bufferEntry = null;
        bufferView = null;

        force.setNodes(nodes);

        if (nodesRef && nodeCount > 0) {
            bufferEntry = retainNodeBuffer(nodesRef, nodeCount);
            if (bufferEntry.ptr) {
                force.setNodeBuffer(bufferEntry.ptr, nodeCount);
                bufferEntry.count = nodeCount;
                bufferView = refreshEntryView(bufferEntry);
            }
        } else {
            force.setNodeBuffer(0, 0);
        }

        return forceWrapper;
    };

    forceWrapper.dispose = function () {
        releaseBuffer();
        nodesRef = null;
        nodeCount = 0;
    };

    forceWrapper.strength = function (_) {
        if (arguments.length) {
            force.setStrength(_);
            return forceWrapper;
        }
        return force.getStrength();
    };

    forceWrapper.distanceMin = function (_) {
        if (arguments.length) {
            force.setDistanceMin(_);
            return forceWrapper;
        }
        return force.getDistanceMin();
    };

    forceWrapper.distanceMax = function (_) {
        if (arguments.length) {
            force.setDistanceMax(_);
            return forceWrapper;
        }
        return force.getDistanceMax();
    };

    forceWrapper.theta = function (_) {
        if (arguments.length) {
            force.setTheta(_);
            return forceWrapper;
        }
        return force.getTheta();
    };

    return forceWrapper;
}


function createForceCollide() {
    ensureModuleReady();

    const force = moduleInstance.createForceCollide();
    let nodesRef = null;
    let nodeCount = 0;
    let bufferEntry = null;
    let bufferView = null;
    let radiusAccessor = null;
    let radiusConstant = 1;
    let radiusPtr = 0;
    let radiusCapacity = 0;
    let radiusView = null;

    const syncNodesToBuffer = (alpha) => {
        if (!nodesRef || nodeCount === 0 || !bufferEntry) return;
        if (bufferEntry.lastUploadAlpha === alpha) return;
        bufferEntry.count = nodeCount;
        bufferView = refreshEntryView(bufferEntry);
        if (!bufferView) return;
        writeNodesToBuffer(nodesRef, bufferView, nodeCount);
        bufferEntry.lastUploadAlpha = alpha;
    };

    const syncBufferToNodes = (_alpha) => {
        if (!nodesRef || nodeCount === 0 || !bufferEntry) return;
        bufferView = refreshEntryView(bufferEntry);
        if (!bufferView) return;
        writeBufferToNodes(bufferView, nodesRef, nodeCount);
    };

    const refreshRadiusView = (length) => {
        if (!radiusPtr || length <= 0) {
            radiusView = null;
            return null;
        }
        if (!radiusView || radiusView.buffer !== moduleInstance.HEAPF32.buffer || radiusView.length !== length) {
            radiusView = new Float32Array(moduleInstance.HEAPF32.buffer, radiusPtr, length);
        }
        return radiusView;
    };

    const ensureRadiusBuffer = (count) => {
        if (count <= 0) {
            return releaseRadiusBuffer();
        }
        if (!radiusPtr || count > radiusCapacity) {
            if (radiusPtr) {
                moduleInstance._free(radiusPtr);
            }
            radiusPtr = moduleInstance._malloc(count * BYTES_PER_VALUE);
            radiusCapacity = count;
        }
        refreshRadiusView(count);
    };

    function releaseRadiusBuffer() {
        if (radiusPtr) {
            moduleInstance._free(radiusPtr);
        }
        radiusPtr = 0;
        radiusCapacity = 0;
        radiusView = null;
        force.setRadiusBuffer(0, 0);
    }

    const populateRadiusBuffer = () => {
        if (!radiusAccessor || !nodesRef || nodeCount === 0) {
            releaseRadiusBuffer();
            return;
        }
        ensureRadiusBuffer(nodeCount);
        if (!radiusView) return;
        for (let i = 0; i < nodeCount; ++i) {
            const node = nodesRef[i];
            const value = radiusAccessor.call(node, node, i, nodesRef);
            radiusView[i] = Number.isFinite(value) && value > 0 ? value : 0;
        }
        force.setRadiusBuffer(radiusPtr, nodeCount);
    };

    function releaseBuffer() {
        if (bufferEntry && nodesRef) {
            force.setNodeBuffer(0, 0);
            releaseNodeBuffer(nodesRef);
        }
        bufferEntry = null;
        bufferView = null;
        releaseRadiusBuffer();
    }

    function forceWrapper(alpha) {
        if (!nodesRef || nodeCount === 0 || !bufferEntry) return;
        syncNodesToBuffer(alpha);
        force.force(alpha ?? 0);
        syncBufferToNodes(alpha);
    }

    forceWrapper.initialize = function (nodes) {
        releaseBuffer();
        nodesRef = nodes || null;
        nodeCount = nodesRef ? nodesRef.length : 0;
        bufferEntry = null;
        bufferView = null;

        force.setNodes(nodes);

        if (nodesRef && nodeCount > 0) {
            bufferEntry = retainNodeBuffer(nodesRef, nodeCount);
            if (bufferEntry.ptr) {
                force.setNodeBuffer(bufferEntry.ptr, nodeCount);
                bufferEntry.count = nodeCount;
                bufferView = refreshEntryView(bufferEntry);
            }
            if (radiusAccessor) {
                populateRadiusBuffer();
            } else {
                force.setRadius(radiusConstant);
                releaseRadiusBuffer();
            }
        } else {
            force.setNodeBuffer(0, 0);
            releaseRadiusBuffer();
        }

        return forceWrapper;
    };

    forceWrapper.dispose = function () {
        releaseBuffer();
        nodesRef = null;
        nodeCount = 0;
    };

    forceWrapper.radius = function (_) {
        if (arguments.length) {
            if (typeof _ === 'function') {
                radiusAccessor = _;
                force.setRadius(0);
                if (nodesRef && nodeCount > 0) {
                    populateRadiusBuffer();
                }
            } else {
                radiusAccessor = null;
                radiusConstant = Number(_) || 0;
                releaseRadiusBuffer();
                force.setRadius(radiusConstant);
            }
            return forceWrapper;
        }
        return radiusAccessor || force.getRadius();
    };

    forceWrapper.strength = function (_) {
        if (arguments.length) {
            force.setStrength(Number(_));
            return forceWrapper;
        }
        return force.getStrength();
    };

    forceWrapper.iterations = function (_) {
        if (arguments.length) {
            const value = Math.max(1, Math.round(Number(_)));
            force.setIterations(value);
            return forceWrapper;
        }
        return force.getIterations();
    };

    return forceWrapper;
}

export const forceX = createAxisForce(x => moduleInstance.createForceX(x), 'x');
export const forceY = createAxisForce(y => moduleInstance.createForceY(y), 'y');
export const forceManyBody = createForceManyBody;
export const forceCollide = createForceCollide;

export function ensureInitialized() {
    return initializationPromise;
}
