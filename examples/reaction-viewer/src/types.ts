export interface Vec3 {
  x: number
  y: number
  z: number
}

export interface AtomData {
  element: number   // atomic number
  position: Vec3
}

export interface BondData {
  i: number
  j: number
  order: number     // 1=single, 2=double, 3=triple
}

export interface MoleculeFrame {
  atoms: AtomData[]
  bonds: BondData[]
  energy: number    // relative energy in kcal/mol or a.u.
}

export type Algorithm = 'uff' | 'pm3' | 'xtb' | 'gfn1' | 'gfn2' | 'hf3c'

export interface ReactionDefinition {
  id: string
  name: string
  description: string
  smirks: string
  reactantSmiles: string
  productSmiles: string
  /** Pre-computed frames for immediate display */
  frames: MoleculeFrame[]
}
