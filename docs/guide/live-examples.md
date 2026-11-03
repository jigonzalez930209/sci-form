---
title: Live examples
---

# Live examples

Every function of the WebAssembly package has a runnable example here — one per function.
Press **Run** and the calculation happens in *your* browser, with the compiled `sci-form-wasm` package built
for this release (the same code as the npm package, in its single-threaded build); nothing is sent to a server.

Each example gets the same starting variables, so the code stays short:

| Variable | Value |
|---|---|
| `sci` | the compiled module |
| `smiles` | `"CCO"` (ethanol) |
| `elements`, `coords` | JSON strings with the atomic numbers and the flat `x, y, z` (Å) of ethanol, embedded with seed 42 |
| `waterEl`, `waterCo` | the same for a fixed water geometry (used by the slower methods) |
| `J` | `JSON.stringify` |

Every function returns a JSON string, so most examples simply `return sci.<function>(…)` and the
page prints the parsed result. Use the same calls in your own code:

```ts
import init, * as sci from 'sci-form-wasm'
await init()
const mol = JSON.parse(sci.embed('CCO', 42))
const pm3 = JSON.parse(sci.compute_pm3(JSON.stringify(mol.elements), JSON.stringify(mol.coords)))
```

<LiveExamples />
