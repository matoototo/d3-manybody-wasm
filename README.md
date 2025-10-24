# d3-manybody-wasm

WebAssembly-powered versions of the most expensive `d3-force` functions:

- `forceManyBody()` (Barnes–Hut charge)
- `forceX()` / `forceY()` (axis forces)
- `forceCollide()` (node collisions)

They’re written in C++ and compiled to WASM. It uses a shared typed array for `[x, y, vx, vy]`, so ticks run inside WASM with minimal JS<->WASM overhead. SIMD is enabled where supported.

## Install

```
npm install d3-manybody-wasm
```

## Usage

```js
import * as d3 from 'd3-force';
import { ensureInitialized, forceManyBody, forceX, forceY, forceCollide } from 'd3-manybody-wasm';

await ensureInitialized();

const sim = d3.forceSimulation(nodes)
  .force('link', d3.forceLink(links).id(d => d.id))
  .force('charge', forceManyBody().strength(-80))
  .force('x', forceX(width / 2))
  .force('y', forceY(height / 2))
  .force('collide', forceCollide().radius(d => d.r));
```

The functions are drop‑in replacements for D3’s forces and follow the same APIs.

## Performance

I observed roughly 3x speedups on medium-large graphs, but your mileage may vary.

## License

MIT — see [LICENSE](LICENSE).
