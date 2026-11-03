// Shared by the documentation site (browser) and scripts/test-live-examples.mjs
// (Node): builds the variables every live example can use and runs one snippet.
//
//   sci       the compiled sci-form WASM module
//   smiles    "CCO" (ethanol)
//   elements  JSON string, atomic numbers of ethanol
//   coords    JSON string, flat x,y,z of ethanol (embedded with seed 42, in Å)
//   waterEl / waterCo   water, a fixed geometry for the slower methods
//   J         JSON.stringify

export const WATER_ELEMENTS = '[8,1,1]';
export const WATER_COORDS = '[0,0,0.1173,0,0.7572,-0.4692,0,-0.7572,-0.4692]';

export function makeScope(sci) {
  const mol = JSON.parse(sci.embed('CCO', 42));
  return {
    sci,
    smiles: 'CCO',
    elements: JSON.stringify(mol.elements),
    coords: JSON.stringify(mol.coords),
    waterEl: WATER_ELEMENTS,
    waterCo: WATER_COORDS,
    J: JSON.stringify,
  };
}

const AsyncFunction = Object.getPrototypeOf(async function () {}).constructor;

/** Run a snippet; returns { value, ms }. A JSON string result is parsed. */
export async function runExample(scope, code) {
  const names = Object.keys(scope);
  const fn = new AsyncFunction(...names, code);
  const t0 = typeof performance !== 'undefined' ? performance.now() : Date.now();
  let value = await fn(...names.map((n) => scope[n]));
  const ms = (typeof performance !== 'undefined' ? performance.now() : Date.now()) - t0;
  if (typeof value === 'string') {
    try {
      value = JSON.parse(value);
    } catch {
      /* plain text result */
    }
  }
  return { value, ms };
}
