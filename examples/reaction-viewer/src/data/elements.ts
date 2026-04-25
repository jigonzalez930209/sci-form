/* ── Element data: CPK colors, covalent radii, names ─────── */

export interface ElementInfo {
  symbol: string
  name: string
  color: string
  radius: number   // covalent radius in Å
  vdwRadius: number // van der Waals radius in Å
}

const ELEMENTS: Record<number, ElementInfo> = {
  1:  { symbol: 'H',  name: 'Hydrogen',  color: '#FFFFFF', radius: 0.31, vdwRadius: 1.20 },
  5:  { symbol: 'B',  name: 'Boron',     color: '#FFB5B5', radius: 0.84, vdwRadius: 1.92 },
  6:  { symbol: 'C',  name: 'Carbon',    color: '#909090', radius: 0.76, vdwRadius: 1.70 },
  7:  { symbol: 'N',  name: 'Nitrogen',  color: '#3050F8', radius: 0.71, vdwRadius: 1.55 },
  8:  { symbol: 'O',  name: 'Oxygen',    color: '#FF0D0D', radius: 0.66, vdwRadius: 1.52 },
  9:  { symbol: 'F',  name: 'Fluorine',  color: '#90E050', radius: 0.57, vdwRadius: 1.47 },
  15: { symbol: 'P',  name: 'Phosphorus',color: '#FF8000', radius: 1.07, vdwRadius: 1.80 },
  16: { symbol: 'S',  name: 'Sulfur',    color: '#FFFF30', radius: 1.05, vdwRadius: 1.80 },
  17: { symbol: 'Cl', name: 'Chlorine',  color: '#1FF01F', radius: 1.02, vdwRadius: 1.75 },
  35: { symbol: 'Br', name: 'Bromine',   color: '#A62929', radius: 1.20, vdwRadius: 1.85 },
  53: { symbol: 'I',  name: 'Iodine',    color: '#940094', radius: 1.39, vdwRadius: 1.98 },
  22: { symbol: 'Ti', name: 'Titanium',  color: '#BFC2C7', radius: 1.60, vdwRadius: 2.00 },
  26: { symbol: 'Fe', name: 'Iron',      color: '#E06633', radius: 1.32, vdwRadius: 2.00 },
  29: { symbol: 'Cu', name: 'Copper',    color: '#C88033', radius: 1.32, vdwRadius: 1.40 },
  30: { symbol: 'Zn', name: 'Zinc',      color: '#7D80B0', radius: 1.22, vdwRadius: 1.39 },
}

export function getElement(z: number): ElementInfo {
  return ELEMENTS[z] ?? {
    symbol: '?',
    name: 'Unknown',
    color: '#FF69B4',
    radius: 1.0,
    vdwRadius: 1.8,
  }
}

/** Bond threshold: two atoms are bonded if distance < (r1 + r2) * factor */
export const BOND_FACTOR = 1.3
