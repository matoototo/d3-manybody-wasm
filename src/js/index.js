import createModule from '../../dist/d3-force-wasm.js';

let moduleInstance = null;
let initializationPromise = null;

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

function createForceWrapper(createForceFunc, coordinateName) {
  return function(coordinate) {
    if (!moduleInstance) {
      throw new Error("WASM module not initialized. Please wait for initialization to complete.");
    }

    const force = createForceFunc(coordinate);

    function forceWrapper(alpha) {
      force.force(alpha);
    }

    forceWrapper.initialize = function(_) {
      force.setNodes(_);
      return forceWrapper;
    };

    forceWrapper.strength = function(_) {
      if (arguments.length) {
        force.setStrength(_);
        return forceWrapper;
      } else {
        return force.getStrength();
      }
    };

    forceWrapper[coordinateName] = function(_) {
      if (arguments.length) {
        force.setCoordinate(_);
        return forceWrapper;
      } else {
        return force.getCoordinate();
      }
    };

    return forceWrapper;
  };
}

export const forceX = createForceWrapper(x => moduleInstance.createForceX(x), 'x');
export const forceY = createForceWrapper(y => moduleInstance.createForceY(y), 'y');

function createForceManyBody() {
  if (!moduleInstance) {
    throw new Error("WASM module not initialized. Please wait for initialization to complete.");
  }

  const force = moduleInstance.createForceManyBody();

  function forceWrapper(alpha) {
    force.force(alpha);
  }

  forceWrapper.initialize = function(nodes) {
    force.setNodes(nodes);
    return forceWrapper;
  };

  forceWrapper.strength = function(_) {
    if (arguments.length) {
      force.setStrength(_);
      return forceWrapper;
    } else {
      return force.getStrength();
    }
  };

  forceWrapper.distanceMin = function(_) {
    if (arguments.length) {
      force.setDistanceMin(_);
      return forceWrapper;
    } else {
      return force.getDistanceMin();
    }
  };

  forceWrapper.distanceMax = function(_) {
    if (arguments.length) {
      force.setDistanceMax(_);
      return forceWrapper;
    } else {
      return force.getDistanceMax();
    }
  };

  forceWrapper.theta = function(_) {
    if (arguments.length) {
      force.setTheta(_);
      return forceWrapper;
    } else {
      return force.getTheta();
    }
  };

  forceWrapper.randomSource = function(_) {
    if (arguments.length) {
      force.setRandom(_);
      return forceWrapper;
    } else {
      return force.getRandom();
    }
  };

  return forceWrapper;
}

export const forceManyBody = createForceManyBody;

export function ensureInitialized() {
  return initializationPromise;
}
