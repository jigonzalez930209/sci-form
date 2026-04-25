import { useRef, useMemo } from 'react'
import { useFrame } from '@react-three/fiber'
import * as THREE from 'three'
import type { AtomData } from '../types'
import { getElement } from '../data/elements'

interface AtomProps {
  atom: AtomData
  index: number
  scale?: number
  selected?: boolean
}

/** Render a single atom as a sphere with CPK coloring */
export default function Atom({ atom, scale = 0.35, selected }: AtomProps) {
  const meshRef = useRef<THREE.Mesh>(null)
  const info = getElement(atom.element)

  const color = useMemo(() => new THREE.Color(info.color), [info.color])
  const r = info.radius * scale

  // Subtle pulse on selection
  useFrame((_, delta) => {
    if (meshRef.current && selected) {
      meshRef.current.scale.setScalar(1 + 0.05 * Math.sin(Date.now() * 0.005))
    }
  })

  return (
    <mesh
      ref={meshRef}
      position={[atom.position.x, atom.position.y, atom.position.z]}
    >
      <sphereGeometry args={[r, 32, 32]} />
      <meshStandardMaterial
        color={color}
        roughness={0.3}
        metalness={0.1}
        emissive={selected ? color : undefined}
        emissiveIntensity={selected ? 0.3 : 0}
      />
    </mesh>
  )
}
