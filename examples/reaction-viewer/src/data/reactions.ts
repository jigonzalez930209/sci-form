import type { ReactionDefinition, MoleculeFrame, AtomData, BondData, Vec3 } from '../types'
import { getElement, BOND_FACTOR } from './elements'

/* ── Helpers ───────────────────────────────────────────────── */

function v(x: number, y: number, z: number): Vec3 { return { x, y, z } }

function lerp3(a: Vec3, b: Vec3, t: number): Vec3 {
  return { x: a.x + (b.x - a.x) * t, y: a.y + (b.y - a.y) * t, z: a.z + (b.z - a.z) * t }
}

function dist(a: Vec3, b: Vec3): number {
  return Math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)
}

/** Detect bonds from elements + positions using covalent radii */
function detectBonds(atoms: AtomData[]): BondData[] {
  const bonds: BondData[] = []
  for (let i = 0; i < atoms.length; i++) {
    for (let j = i + 1; j < atoms.length; j++) {
      const ri = getElement(atoms[i].element).radius
      const rj = getElement(atoms[j].element).radius
      const d = dist(atoms[i].position, atoms[j].position)
      if (d < (ri + rj) * BOND_FACTOR && d > 0.4) {
        bonds.push({ i, j, order: 1 })
      }
    }
  }
  return bonds
}

/** Generate interpolated frames with energy profile */
function generateFrames(
  elements: number[],
  startPos: Vec3[],
  endPos: Vec3[],
  nFrames: number,
  energyFn: (t: number) => number,
): MoleculeFrame[] {
  const frames: MoleculeFrame[] = []
  for (let f = 0; f < nFrames; f++) {
    const t = f / (nFrames - 1)
    const atoms: AtomData[] = elements.map((el, i) => ({
      element: el,
      position: lerp3(startPos[i], endPos[i], t),
    }))
    const bonds = detectBonds(atoms)
    frames.push({ atoms, bonds, energy: energyFn(t) })
  }
  return frames
}

/* ── Reaction 1: H₂ Dissociation ──────────────────────────── */
function makeH2Dissociation(nFrames: number): MoleculeFrame[] {
  const elements = [1, 1]
  const start = [v(-0.37, 0, 0), v(0.37, 0, 0)]
  const end = [v(-1.8, 0, 0), v(1.8, 0, 0)]
  // Symmetric double-well: V(t) = 4t²(1-t)² * barrier, shifted to peak at t=0.5
  return generateFrames(elements, start, end, nFrames, t => {
    const s = 2 * t - 1
    return (1 - s * s) * 110 // barrier ~110 kcal/mol
  })
}

/* ── Reaction 2: SN2 — CH₃F + Cl⁻ → CH₃Cl + F⁻ ─────────── */
function makeSN2(nFrames: number): MoleculeFrame[] {
  // atoms: C(0), H(1), H(2), H(3), F(4), Cl(5)
  const elements = [6, 1, 1, 1, 9, 17]
  const start: Vec3[] = [
    v(0, 0, 0),         // C
    v(0, 1.09, 0),      // H up
    v(1.03, -0.36, 0),  // H right-front
    v(-0.51, -0.36, 0.89), // H left-back
    v(-0.51, -0.36, -0.89), // F (leaving group, close)
    v(3.5, 0, 0),       // Cl⁻ incoming (far)
  ]
  // Walden inversion: C inverts, Cl bonds, F leaves
  const end: Vec3[] = [
    v(0, 0, 0),         // C
    v(0, -1.09, 0),     // H flipped
    v(-1.03, 0.36, 0),  // H flipped
    v(0.51, 0.36, 0.89),  // H flipped
    v(-3.5, 0, 0),      // F⁻ leaves (far)
    v(0.51, 0.36, -0.89), // Cl bonded (close) — actually let me fix
  ]
  // Corrections: Cl comes in, F goes out, H's invert
  const endFixed: Vec3[] = [
    v(0, 0, 0),
    v(0, -1.09, 0),
    v(-1.03, 0.36, 0),
    v(0.51, 0.36, 0.89),
    v(-3.8, 0, 0),      // F⁻ far left
    v(1.78, 0, 0),      // Cl bonded right
  ]
  return generateFrames(elements, start, endFixed, nFrames, t => {
    // Asymmetric barrier peaking around t=0.45
    const peak = 0.45
    const w = 0.25
    const barrier = 15 // kcal/mol
    const rxnEnergy = -8 // exothermic
    return barrier * Math.exp(-(((t - peak) / w) ** 2)) + rxnEnergy * t
  })
}

/* ── Reaction 3: Carboxylic acid deprotonation ────────────── */
function makeDeprotonation(nFrames: number): MoleculeFrame[] {
  // Acetic acid: CH₃COOH → CH₃COO⁻ + H⁺
  // atoms: C(0) C(1) O(2) O(3) H(4) H(5) H(6) H(7-leaving)
  const elements = [6, 6, 8, 8, 1, 1, 1, 1]
  const start: Vec3[] = [
    v(0, 0, 0),          // C methyl
    v(1.52, 0, 0),       // C carbonyl
    v(2.15, 1.08, 0),    // O=
    v(2.15, -1.08, 0),   // O-H
    v(-0.36, 0, 1.03),   // H
    v(-0.36, 1.03, -0.36), // H
    v(-0.36, -1.03, -0.36), // H
    v(3.10, -1.08, 0),   // H (acidic, leaving)
  ]
  const end: Vec3[] = [
    v(0, 0, 0),
    v(1.52, 0, 0),
    v(2.15, 1.08, 0),    // O⁻ (symmetric)
    v(2.15, -1.08, 0),   // O⁻ (symmetric)
    v(-0.36, 0, 1.03),
    v(-0.36, 1.03, -0.36),
    v(-0.36, -1.03, -0.36),
    v(5.0, -1.5, 0),     // H⁺ departs
  ]
  return generateFrames(elements, start, end, nFrames, t => {
    // Endothermic: energy rises monotonically with small barrier
    return 12 * t + 5 * Math.exp(-(((t - 0.6) / 0.2) ** 2))
  })
}

/* ── Reaction 4: Ketone reduction (acetone + H) ──────────── */
function makeKetoneReduction(nFrames: number): MoleculeFrame[] {
  // Simplified: C=O + H → C-OH
  // Acetone: (CH₃)₂C=O + H → (CH₃)₂CHOH
  // Atoms: C(0-central), C(1-methyl), C(2-methyl), O(3), H(4-6 methyl), H(7-9 methyl), H(10-incoming)
  const elements = [6, 6, 6, 8, 1, 1, 1, 1, 1, 1, 1]
  const start: Vec3[] = [
    v(0, 0, 0),           // C central (sp2)
    v(-1.27, 0.73, 0),    // CH₃ left
    v(1.27, 0.73, 0),     // CH₃ right
    v(0, -1.21, 0),       // O (double bond)
    v(-1.8, 0.3, 0.9),    // H
    v(-1.8, 0.3, -0.9),   // H
    v(-1.27, 1.82, 0),    // H
    v(1.8, 0.3, 0.9),     // H
    v(1.8, 0.3, -0.9),    // H
    v(1.27, 1.82, 0),     // H
    v(0, -3.5, 0),        // H incoming (far)
  ]
  const end: Vec3[] = [
    v(0, 0, 0),           // C central (sp3 now)
    v(-1.27, 0.73, 0),
    v(1.27, 0.73, 0),
    v(0, -1.43, 0),       // O (single bond, longer)
    v(-1.8, 0.3, 0.9),
    v(-1.8, 0.3, -0.9),
    v(-1.27, 1.82, 0),
    v(1.8, 0.3, 0.9),
    v(1.8, 0.3, -0.9),
    v(1.27, 1.82, 0),
    v(0, -2.38, 0),       // H bonded to O
  ]
  return generateFrames(elements, start, end, nFrames, t => {
    // Exothermic with moderate barrier
    const barrier = 25
    const rxn = -15
    return barrier * 4 * t * (1 - t) + rxn * t
  })
}

/* ── Reaction 5: Diels-Alder (butadiene + ethylene) ──────── */
function makeDielsAlder(nFrames: number): MoleculeFrame[] {
  // Butadiene (C₁=C₂-C₃=C₄) + Ethylene (C₅=C₆) → Cyclohexene
  // 6 C + 8 H = 14 atoms, but simplified to just the carbons + key H's
  // C(0..5), H(6..13)
  const elements = [6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1]
  // Butadiene (s-cis conformation) in xy-plane
  const start: Vec3[] = [
    v(-1.8, 1.0, 0),    // C1
    v(-0.6, 0.5, 0),    // C2
    v(0.6, 0.5, 0),     // C3
    v(1.8, 1.0, 0),     // C4
    v(-0.6, -2.8, 0),   // C5 (ethylene, far below)
    v(0.6, -2.8, 0),    // C6
    // Hydrogens on butadiene
    v(-2.7, 0.5, 0),    // H on C1
    v(-1.8, 2.1, 0),    // H on C1
    v(-0.6, -0.6, 0),   // H on C2 — no, C2 is internal... let me fix
    v(0.6, -0.6, 0),    // H on C3
    v(1.8, 2.1, 0),     // H on C4
    v(2.7, 0.5, 0),     // H on C4
    // Hydrogens on ethylene
    v(-0.6, -3.9, 0),   // H on C5
    v(0.6, -3.9, 0),    // H on C6
  ]
  // Cyclohexene (C1=C2-C3-C6-C5-C4, ring closed)
  // Chair-like, C1-C4 bond and C5-C6 stay bonded
  const end: Vec3[] = [
    v(-1.2, 0.7, -0.3), // C1  
    v(0, 1.4, 0),       // C2
    v(1.2, 0.7, 0.3),   // C3
    v(1.2, -0.7, -0.3), // C4 (now bonded to C5 neighbor)
    v(0, -1.4, 0),      // C5
    v(-1.2, -0.7, 0.3), // C6 (now bonded to C1 neighbor)
    // H's move with their carbons
    v(-2.1, 1.2, -0.5),
    v(-1.2, 0.7, -1.4),
    v(0, 2.5, 0),
    v(2.1, 1.2, 0.5),
    v(1.2, -0.7, -1.4),
    v(2.1, -1.2, -0.5),
    v(0, -2.5, 0),
    v(-2.1, -1.2, 0.5),
  ]
  return generateFrames(elements, start, end, nFrames, t => {
    // Concerted, moderate barrier, exothermic
    const barrier = 30
    const rxn = -40
    return barrier * 4 * t * (1 - t) + rxn * t
  })
}

/* ── Export all 5 reactions ────────────────────────────────── */

export function getReactions(nFrames: number = 21): ReactionDefinition[] {
  return [
    {
      id: 'h2-dissociation',
      name: 'H₂ Dissociation',
      description: 'Homolytic dissociation of molecular hydrogen. Symmetric double-well potential with ΔE‡ ≈ 110 kcal/mol.',
      smirks: '[H:1][H:2]>>[H:1].[H:2]',
      reactantSmiles: '[H][H]',
      productSmiles: '[H].[H]',
      frames: makeH2Dissociation(nFrames),
    },
    {
      id: 'sn2',
      name: 'SN2 Halide Exchange',
      description: 'Bimolecular nucleophilic substitution: Cl⁻ attacks CH₃F with backside approach and Walden inversion.',
      smirks: '[C:1][F:2].[Cl-:3]>>[C:1][Cl:3].[F-:2]',
      reactantSmiles: 'CF.[Cl-]',
      productSmiles: 'CCl.[F-]',
      frames: makeSN2(nFrames),
    },
    {
      id: 'deprotonation',
      name: 'Acid Deprotonation',
      description: 'Carboxylic acid deprotonation: CH₃COOH → CH₃COO⁻ + H⁺. Endothermic proton departure.',
      smirks: '[C:1](=O)[OH:2]>>[C:1](=O)[O-:2]',
      reactantSmiles: 'CC(=O)O',
      productSmiles: 'CC(=O)[O-]',
      frames: makeDeprotonation(nFrames),
    },
    {
      id: 'ketone-reduction',
      name: 'Ketone Reduction',
      description: 'Simplified reduction of acetone: C=O + H → C-OH. Exothermic with moderate barrier.',
      smirks: '[C:1]=[O:2]>>[C:1][OH:2]',
      reactantSmiles: 'CC(=O)C',
      productSmiles: 'CC(O)C',
      frames: makeKetoneReduction(nFrames),
    },
    {
      id: 'diels-alder',
      name: 'Diels-Alder',
      description: 'Concerted [4+2] cycloaddition: butadiene + ethylene → cyclohexene. Suprafacial approach.',
      smirks: '[C:1]=[C:2][C:3]=[C:4].[C:5]=[C:6]>>[C:1]1[C:2]=[C:3][C:4][C:6][C:5]1',
      reactantSmiles: 'C=CC=C.C=C',
      productSmiles: 'C1CC=CCC1',
      frames: makeDielsAlder(nFrames),
    },
  ]
}
