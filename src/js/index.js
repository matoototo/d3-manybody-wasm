import createModule from '../../dist/d3-manybody-wasm.js';

let moduleInstance = null;
let initializationPromise = null;

const STRIDE = 4;
const BYTES_PER_VALUE = Float64Array.BYTES_PER_ELEMENT;
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
    entry = { ptr: 0, capacity: 0, view: null, refCount: 0, count: 0 };
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

    const syncNodesToBuffer = () => {
      if (!nodesRef || nodeCount === 0 || !bufferEntry) return;
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

    function forceWrapper(alpha) {
      if (!nodesRef || nodeCount === 0 || !bufferEntry) return;
      syncNodesToBuffer();
      force.force(alpha);
      syncBufferToNodes();
    }

    forceWrapper.initialize = function(nodes) {
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
          syncNodesToBuffer();
        }
      }

      if (!bufferEntry) {
        force.setNodeBuffer(0, 0);
      }

      return forceWrapper;
    };

    forceWrapper.dispose = function() {
      releaseBuffer();
      nodesRef = null;
      nodeCount = 0;
    };

    forceWrapper.strength = function(_) {
      if (arguments.length) {
        force.setStrength(_);
        return forceWrapper;
      }
      return force.getStrength();
    };

    forceWrapper[coordinateName] = function(_) {
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

  const syncNodesToBuffer = () => {
    if (!nodesRef || nodeCount === 0 || !bufferEntry) return;
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

  function forceWrapper(alpha) {
    if (!nodesRef || nodeCount === 0 || !bufferEntry) return;
    syncNodesToBuffer();
    force.force(alpha);
    syncBufferToNodes();
  }

  forceWrapper.initialize = function(nodes) {
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
        syncNodesToBuffer();
      }
    } else {
      force.setNodeBuffer(0, 0);
    }

    return forceWrapper;
  };

  forceWrapper.dispose = function() {
    releaseBuffer();
    nodesRef = null;
    nodeCount = 0;
  };

  forceWrapper.strength = function(_) {
    if (arguments.length) {
      force.setStrength(_);
      return forceWrapper;
    }
    return force.getStrength();
  };

  forceWrapper.distanceMin = function(_) {
    if (arguments.length) {
      force.setDistanceMin(_);
      return forceWrapper;
    }
    return force.getDistanceMin();
  };

  forceWrapper.distanceMax = function(_) {
    if (arguments.length) {
      force.setDistanceMax(_);
      return forceWrapper;
    }
    return force.getDistanceMax();
  };

  forceWrapper.theta = function(_) {
    if (arguments.length) {
      force.setTheta(_);
      return forceWrapper;
    }
    return force.getTheta();
  };

  return forceWrapper;
}

export const forceX = createAxisForce(x => moduleInstance.createForceX(x), 'x');
export const forceY = createAxisForce(y => moduleInstance.createForceY(y), 'y');
export const forceManyBody = createForceManyBody;

export function ensureInitialized() {
  return initializationPromise;
}
