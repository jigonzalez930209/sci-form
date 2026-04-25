import type { MoleculeFrame } from '../types'
import Atom from './Atom'
import Bond from './Bond'

interface MoleculeProps {
  frame: MoleculeFrame
  atomScale?: number
  bondRadius?: number
}

/** Render a complete molecule for one frame: atoms + bonds */
export default function Molecule({ frame, atomScale = 0.35, bondRadius = 0.08 }: MoleculeProps) {
  return (
    <group>
      {frame.bonds.map((bond, idx) => (
        <Bond key={`b-${idx}`} bond={bond} atoms={frame.atoms} radius={bondRadius} />
      ))}
      {frame.atoms.map((atom, idx) => (
        <Atom key={`a-${idx}`} atom={atom} index={idx} scale={atomScale} />
      ))}
    </group>
  )
}
