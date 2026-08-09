import createModule from '../../dist/d3-manybody-wasm.js';

let moduleInstance = null;
let initializationPromise = null;

const STRIDE = 4;
const NODE_BYTES_PER_VALUE = Float64Array.BYTES_PER_ELEMENT;
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
    const heapBuffer = moduleInstance.HEAPF64.buffer;

    if (!entry.view || entry.view.buffer !== heapBuffer || entry.view.length !== desiredLength) {
        entry.view = new Float64Array(heapBuffer, entry.ptr, desiredLength);
    }

    return entry.view;
}

function retainNodeBuffer(nodes, count) {
    ensureModuleReady();

    let entry = nodeBufferRegistry.get(nodes);
    if (!entry) {
        entry = { ptr: 0, capacity: 0, view: null, refCount: 0, count: 0, syncAlpha: Number.NaN, seenForces: new Set() };
        nodeBufferRegistry.set(nodes, entry);
    }

    const required = count * STRIDE;
    if (required > entry.capacity || entry.ptr === 0) {
        if (entry.ptr) {
            moduleInstance._free(entry.ptr);
        }
        if (required > 0) {
            entry.ptr = moduleInstance._malloc(required * NODE_BYTES_PER_VALUE);
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
    entry.syncAlpha = Number.NaN;
    entry.seenForces.clear();
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
        const node = nodes[i];
        const base = i * STRIDE;
        view[base] = node.x;
        view[base + 1] = node.y;
        view[base + 2] = node.vx;
        view[base + 3] = node.vy;
    }
}

function writeBufferToNodes(view, nodes, count) {
    if (!view) return;
    for (let i = 0; i < count; ++i) {
        const node = nodes[i];
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

function createAxisForce(createForceFunc, coordinateName) {
    return function axisForceFactory(coordinate) {
        ensureModuleReady();

        const force = createForceFunc(coordinate);
        let nodesRef = null;
        let nodeCount = 0;
        let bufferEntry = null;
        let bufferView = null;
        const syncToken = {};

        const syncNodesToBuffer = (alpha) => {
            if (!nodesRef || nodeCount === 0 || !bufferEntry) return;
            if (!Object.is(bufferEntry.syncAlpha, alpha)) {
                bufferEntry.syncAlpha = alpha;
                bufferEntry.seenForces.clear();
            } else if (bufferEntry.seenForces.has(syncToken)) {
                // The same force appearing again marks a new tick even when a
                // simulation deliberately keeps alpha constant.
                bufferEntry.seenForces.clear();
            }
            const shouldUpload = bufferEntry.seenForces.size === 0;
            bufferEntry.seenForces.add(syncToken);
            if (!shouldUpload) return;
            bufferEntry.count = nodeCount;
            bufferView = refreshEntryView(bufferEntry);
            if (!bufferView) return;
            writeNodesToBuffer(nodesRef, bufferView, nodeCount);
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
            forceWrapper._prepare(alpha);
            forceWrapper._apply(alpha);
            forceWrapper._flush(alpha);
        }

        forceWrapper._prepare = syncNodesToBuffer;
        forceWrapper._apply = (alpha) => force.force(alpha);
        forceWrapper._flush = syncBufferToNodes;

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

function createForceManyBody(createNativeForce) {
    ensureModuleReady();

    const force = createNativeForce();
    let nodesRef = null;
    let nodeCount = 0;
    let bufferEntry = null;
    let bufferView = null;
    let strengthSetting = -30;
    const syncToken = {};

    const configureStrengths = () => {
        if (typeof strengthSetting !== 'function') {
            force.setStrength(Number(strengthSetting));
            return;
        }
        if (!nodesRef) return;
        const values = new Float32Array(nodeCount);
        for (let index = 0; index < nodeCount; ++index) {
            values[index] = Number(strengthSetting.call(nodesRef[index], nodesRef[index], index, nodesRef));
        }
        force.setStrengths(values);
    };

    const syncNodesToBuffer = (alpha) => {
        if (!nodesRef || nodeCount === 0 || !bufferEntry) return;
        if (!Object.is(bufferEntry.syncAlpha, alpha)) {
            bufferEntry.syncAlpha = alpha;
            bufferEntry.seenForces.clear();
        } else if (bufferEntry.seenForces.has(syncToken)) {
            bufferEntry.seenForces.clear();
        }
        const shouldUpload = bufferEntry.seenForces.size === 0;
        bufferEntry.seenForces.add(syncToken);
        if (!shouldUpload) return;
        bufferEntry.count = nodeCount;
        bufferView = refreshEntryView(bufferEntry);
        if (!bufferView) return;
        writeNodesToBuffer(nodesRef, bufferView, nodeCount);
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
        forceWrapper._prepare(alpha);
        forceWrapper._apply(alpha);
        forceWrapper._flush(alpha);
    }

    forceWrapper._prepare = syncNodesToBuffer;
    forceWrapper._apply = (alpha) => force.force(alpha);
    forceWrapper._flush = syncBufferToNodes;
    if (typeof force.beginPrecompute === 'function') {
        forceWrapper._before = alpha => force.beginPrecompute(alpha);
    }

    forceWrapper.initialize = function (nodes) {
        releaseBuffer();
        nodesRef = nodes || null;
        nodeCount = nodesRef ? nodesRef.length : 0;
        bufferEntry = null;
        bufferView = null;

        force.setNodes(nodes);
        configureStrengths();
        force.prepareWorkers?.();

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
        force.shutdownWorkers?.();
        releaseBuffer();
        nodesRef = null;
        nodeCount = 0;
    };

    forceWrapper.strength = function (_) {
        if (arguments.length) {
            strengthSetting = _;
            configureStrengths();
            return forceWrapper;
        }
        return typeof strengthSetting === 'function' ? strengthSetting : () => Number(strengthSetting);
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

    forceWrapper.axes = function (x, y, strength) {
        if (typeof force.setAxes !== 'function') throw new Error('axes fusion is only available for the legacy force');
        force.setAxes(Number(x), Number(y), Number(strength));
        return forceWrapper;
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
    const syncToken = {};

    const syncNodesToBuffer = (alpha) => {
        if (!nodesRef || nodeCount === 0 || !bufferEntry) return;
        if (!Object.is(bufferEntry.syncAlpha, alpha)) {
            bufferEntry.syncAlpha = alpha;
            bufferEntry.seenForces.clear();
        } else if (bufferEntry.seenForces.has(syncToken)) {
            bufferEntry.seenForces.clear();
        }
        const shouldUpload = bufferEntry.seenForces.size === 0;
        bufferEntry.seenForces.add(syncToken);
        if (!shouldUpload) return;
        bufferEntry.count = nodeCount;
        bufferView = refreshEntryView(bufferEntry);
        if (!bufferView) return;
        writeNodesToBuffer(nodesRef, bufferView, nodeCount);
    };

    const syncBufferToNodes = (_alpha) => {
        if (!nodesRef || nodeCount === 0 || !bufferEntry) return;
        bufferView = refreshEntryView(bufferEntry);
        if (!bufferView) return;
        writeBufferToNodes(bufferView, nodesRef, nodeCount);
    };

    function releaseRadiusBuffer() {
        force.setRadiusBuffer(0, 0);
    }

    const populateRadiusBuffer = () => {
        if (!radiusAccessor || !nodesRef || nodeCount === 0) {
            releaseRadiusBuffer();
            return;
        }
        const values = new Float32Array(nodeCount);
        for (let i = 0; i < nodeCount; ++i) {
            const node = nodesRef[i];
            const value = radiusAccessor.call(node, node, i, nodesRef);
            values[i] = Number(value);
        }
        releaseRadiusBuffer();
        force.setRadii(values);
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
        forceWrapper._prepare(alpha);
        forceWrapper._apply(alpha);
        forceWrapper._flush(alpha);
    }

    forceWrapper._prepare = syncNodesToBuffer;
    forceWrapper._apply = (alpha) => force.force(alpha ?? 0);
    forceWrapper._flush = syncBufferToNodes;

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
                if (nodesRef && nodeCount > 0) populateRadiusBuffer();
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

function createForceLink(initialLinks = []) {
    ensureModuleReady();

    const force = moduleInstance.createForceLink();
    let linksRef = initialLinks || [];
    let nodesRef = null;
    let nodeCount = 0;
    let bufferEntry = null;
    let bufferView = null;
    let idAccessor = (node) => node.index;
    let strengthAccessor = null;
    let distanceAccessor = () => 30;
    let counts = [];
    let iterationCount = 1;
    const syncToken = {};

    const syncNodesToBuffer = (alpha) => {
        if (!nodesRef || nodeCount === 0 || !bufferEntry) return;
        if (!Object.is(bufferEntry.syncAlpha, alpha)) {
            bufferEntry.syncAlpha = alpha;
            bufferEntry.seenForces.clear();
        } else if (bufferEntry.seenForces.has(syncToken)) {
            bufferEntry.seenForces.clear();
        }
        const shouldUpload = bufferEntry.seenForces.size === 0;
        bufferEntry.seenForces.add(syncToken);
        if (!shouldUpload) return;
        bufferEntry.count = nodeCount;
        bufferView = refreshEntryView(bufferEntry);
        if (!bufferView) return;
        writeNodesToBuffer(nodesRef, bufferView, nodeCount);
    };

    const syncBufferToNodes = () => {
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

    function defaultStrength(link) {
        return 1 / Math.min(counts[link.source.index], counts[link.target.index]);
    }

    function configureLinks() {
        if (!nodesRef) return;
        const nodeById = new Map(nodesRef.map((node, index) => [idAccessor(node, index, nodesRef), node]));
        counts = new Array(nodeCount).fill(0);

        for (let index = 0; index < linksRef.length; ++index) {
            const link = linksRef[index];
            link.index = index;
            if (typeof link.source !== 'object') link.source = nodeById.get(link.source);
            if (typeof link.target !== 'object') link.target = nodeById.get(link.target);
            if (!link.source) throw new Error(`node not found: ${link.source}`);
            if (!link.target) throw new Error(`node not found: ${link.target}`);
            counts[link.source.index] += 1;
            counts[link.target.index] += 1;
        }

        const sources = new Int32Array(linksRef.length);
        const targets = new Int32Array(linksRef.length);
        const biases = new Float64Array(linksRef.length);
        const strengths = new Float64Array(linksRef.length);
        const distances = new Float64Array(linksRef.length);
        const getStrength = strengthAccessor || defaultStrength;
        for (let index = 0; index < linksRef.length; ++index) {
            const link = linksRef[index];
            const sourceCount = counts[link.source.index];
            const targetCount = counts[link.target.index];
            sources[index] = link.source.index;
            targets[index] = link.target.index;
            biases[index] = sourceCount / (sourceCount + targetCount);
            strengths[index] = Number(getStrength(link, index, linksRef));
            distances[index] = Number(distanceAccessor(link, index, linksRef));
        }
        force.setLinks(sources, targets, biases, strengths, distances);
    }

    function forceWrapper(alpha) {
        if (!nodesRef || nodeCount === 0 || !bufferEntry) return;
        forceWrapper._prepare(alpha);
        forceWrapper._apply(alpha);
        forceWrapper._flush(alpha);
    }

    forceWrapper._prepare = syncNodesToBuffer;
    forceWrapper._apply = (alpha) => force.force(alpha ?? 0);
    forceWrapper._flush = syncBufferToNodes;

    forceWrapper.initialize = function (nodes) {
        releaseBuffer();
        nodesRef = nodes || null;
        nodeCount = nodesRef ? nodesRef.length : 0;
        if (nodesRef && nodeCount > 0) {
            bufferEntry = retainNodeBuffer(nodesRef, nodeCount);
            force.setNodeBuffer(bufferEntry.ptr, nodeCount);
            bufferView = refreshEntryView(bufferEntry);
            configureLinks();
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

    forceWrapper.links = function (_) {
        if (!arguments.length) return linksRef;
        linksRef = _ || [];
        configureLinks();
        return forceWrapper;
    };

    forceWrapper.id = function (_) {
        if (!arguments.length) return idAccessor;
        idAccessor = _;
        configureLinks();
        return forceWrapper;
    };

    forceWrapper.iterations = function (_) {
        if (!arguments.length) return iterationCount;
        iterationCount = Math.max(1, Math.round(Number(_)));
        force.setIterations(iterationCount);
        return forceWrapper;
    };

    forceWrapper.strength = function (_) {
        if (!arguments.length) return strengthAccessor || defaultStrength;
        strengthAccessor = typeof _ === 'function' ? _ : () => Number(_);
        configureLinks();
        return forceWrapper;
    };

    forceWrapper.distance = function (_) {
        if (!arguments.length) return distanceAccessor;
        distanceAccessor = typeof _ === 'function' ? _ : () => Number(_);
        configureLinks();
        return forceWrapper;
    };

    return forceWrapper;
}

function createForceBundle(forces) {
    const bundledForces = (forces || []).filter(Boolean);
    if (!bundledForces.length || bundledForces.some(force =>
        typeof force._prepare !== 'function' ||
        typeof force._apply !== 'function' ||
        typeof force._flush !== 'function'
    )) {
        throw new TypeError('forceBundle expects initialized d3-manybody-wasm forces');
    }

    const linkForce = bundledForces.find(force => typeof force.links === 'function');
    function bundle(alpha) {
        bundledForces[0]._prepare(alpha);
        for (const force of bundledForces) force._before?.(alpha);
        for (const force of bundledForces) force._apply(alpha);
        bundledForces[bundledForces.length - 1]._flush(alpha);
    }

    bundle.initialize = function (nodes, random) {
        for (const force of bundledForces) force.initialize(nodes, random);
        return bundle;
    };

    bundle.dispose = function () {
        for (const force of bundledForces) force.dispose?.();
    };

    for (const method of ['links', 'id', 'distance', 'iterations']) {
        if (!linkForce || typeof linkForce[method] !== 'function') continue;
        bundle[method] = function (_) {
            if (!arguments.length) return linkForce[method]();
            linkForce[method](_);
            return bundle;
        };
    }

    bundle.forces = () => [...bundledForces];
    return bundle;
}

export const forceX = createAxisForce(x => moduleInstance.createForceX(x), 'x');
export const forceY = createAxisForce(y => moduleInstance.createForceY(y), 'y');
export const forceManyBody = () => createForceManyBody(() => moduleInstance.createForceManyBodyLegacy());
export const forceManyBodyLegacy = forceManyBody;
export const forceCollide = createForceCollide;
export const forceLink = createForceLink;
export const forceBundle = createForceBundle;

export function ensureInitialized() {
    return initializeWasm();
}
