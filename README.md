# d3-manybody-wasm

This module provides WebAssembly-powered implementations of force algorithms used in d3-force: `forceManyBody()`, `forceX()`, and `forceY()`. By leveraging WebAssembly, these implementations can significantly improve the performance of force-directed graph layouts, especially for large graphs, without compromising on layout quality.

This is done through a C++ implementation, which is then compiled to WebAssembly. This allows for faster execution of the force calculations, particularly on larger datasets where JavaScript performance may become a bottleneck.

## Installing

If you use NPM, `npm install d3-manybody-wasm`. Otherwise, clone [latest version](https://github.com/matoototo/d3-manybody-wasm).

## Usage

Before using any of the force functions, you must ensure that the WebAssembly module is initialized. Use the `ensureInitialized()` function to wait for the initialization to complete:

```javascript
import { ensureInitialized, forceManyBody, forceX, forceY } from 'd3-manybody-wasm';

ensureInitialized().then(() => {
  const simulation = d3.forceSimulation(nodes)
    .force("link", d3.forceLink().id(d => d.id))
    .force("charge", forceManyBody())
    .force("x", forceX())
    .force("y", forceY());

  // ... rest of your simulation code
});
```

Aside from that, The `forceManyBody()`, `forceX()`, and `forceY()` functions are designed as drop-in replacements for their d3-force counterparts, using the same API.

## API Reference

### Many-Body Force

<a name="forceManyBody" href="#forceManyBody">#</a> d3.<b>forceManyBody</b>()

Creates a new many-body force with the default parameters. This force can be used to simulate gravity (attraction) or electrostatic charge (repulsion).

<a name="manyBody_strength" href="#manyBody_strength">#</a> <i>manyBody</i>.<b>strength</b>([<i>strength</i>])

If *strength* is specified, sets the strength accessor to the specified number or function and returns this force. If *strength* is not specified, returns the current strength accessor.

<a name="manyBody_theta" href="#manyBody_theta">#</a> <i>manyBody</i>.<b>theta</b>([<i>theta</i>])

If *theta* is specified, sets the Barnes–Hut approximation criterion to the specified number and returns this force. If *theta* is not specified, returns the current value.

<a name="manyBody_distanceMin" href="#manyBody_distanceMin">#</a> <i>manyBody</i>.<b>distanceMin</b>([<i>distance</i>])

If *distance* is specified, sets the minimum distance between nodes and returns this force. If *distance* is not specified, returns the current minimum distance.

<a name="manyBody_distanceMax" href="#manyBody_distanceMax">#</a> <i>manyBody</i>.<b>distanceMax</b>([<i>distance</i>])

If *distance* is specified, sets the maximum distance between nodes and returns this force. If *distance* is not specified, returns the current maximum distance.

### Position Forces

<a name="forceX" href="#forceX">#</a> d3.<b>forceX</b>([<i>x</i>])

Creates a new position force along the x-axis towards the given position *x*.

<a name="forceY" href="#forceY">#</a> d3.<b>forceY</b>([<i>y</i>])

Creates a new position force along the y-axis towards the given position *y*.

Both position forces share the following methods:

<a name="force_strength" href="#force_strength">#</a> <i>force</i>.<b>strength</b>([<i>strength</i>])

If *strength* is specified, sets the strength accessor to the specified number or function and returns this force. If *strength* is not specified, returns the current strength accessor.

<a name="force_x" href="#force_x">#</a> <i>force</i>.<b>x</b>([<i>x</i>])

If *x* is specified, sets the x-coordinate accessor to the specified number or function and returns this force. If *x* is not specified, returns the current x-coordinate accessor.

<a name="force_y" href="#force_y">#</a> <i>force</i>.<b>y</b>([<i>y</i>])

If *y* is specified, sets the y-coordinate accessor to the specified number or function and returns this force. If *y* is not specified, returns the current y-coordinate accessor.

## Performance

The WebAssembly implementation can provide decent performance improvements, especially for larger graphs. However, the exact performance gain may vary depending on the specific use case and hardware. I encourage you to benchmark this implementation against the standard D3 force implementation and see for yourself – in particular, data copying between C++ and JavaScript can be a bottleneck. Setting up shared memory and hacking the package a bit to accommodate it might be a good idea if you're looking to optimize further.


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.