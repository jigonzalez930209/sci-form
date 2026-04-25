import { useMemo } from 'react'
import * as THREE from 'three'
import type { AtomData, BondData } from '../types'

interface BondProps {
  bond: BondData
  atoms: AtomData[]
  radius?: number
}

/** Render a bond as a cylinder between two atoms */
export default function Bond({ bond, atoms, radius = 0.08 }: BondProps) {
  const a = atoms[bond.i].position
  const b = atoms[bond.j].position

  const { position, quaternion, length } = useMemo(() => {
    const start = new THREE.Vector3(a.x, a.y, a.z)
    const end = new THREE.Vector3(b.x, b.y, b.z)
    const mid = new THREE.Vector3().addVectors(start, end).multiplyScalar(0.5)
    const direction = new THREE.Vector3().subVectors(end, start)
    const len = direction.length()
    direction.normalize()

    // Cylinder default axis is Y; rotate to match bond direction
    const quat = new THREE.Quaternion()
    quat.setFromUnitVectors(new THREE.Vector3(0, 1, 0), direction)

    return { position: mid, quaternion: quat, length: len }
  }, [a.x, a.y, a.z, b.x, b.y, b.z])

  if (bond.order === 1) {
    return (
      <mesh position={position} quaternion={quaternion}>
        <cylinderGeometry args={[radius, radius, length, 8]} />
        <meshStandardMaterial color="#888888" roughness={0.4} metalness={0.05} />
      </mesh>
    )
  }

  // Double/triple bonds: offset cylinders
  const offsets = bond.order === 2 ? [-0.06, 0.06] : [-0.08, 0, 0.08]
  const dir = new THREE.Vector3(a.x, a.y, a.z).sub(new THREE.Vector3(b.x, b.y, b.z)).normalize()
  const perp = new THREE.Vector3()
  if (Math.abs(dir.x) < 0.9) perp.crossVectors(dir, new THREE.Vector3(1, 0, 0)).normalize()
  else perp.crossVectors(dir, new THREE.Vector3(0, 1, 0)).normalize()

  return (
    <group>
      {offsets.map((off, idx) => {
        const offset = perp.clone().multiplyScalar(off)
        return (
          <mesh
            key={idx}
            position={[position.x + offset.x, position.y + offset.y, position.z + offset.z]}
            quaternion={quaternion}
          >
            <cylinderGeometry args={[radius * 0.7, radius * 0.7, length, 8]} />
            <meshStandardMaterial color="#888888" roughness={0.4} metalness={0.05} />
          </mesh>
        )
      })}
    </group>
  )
}
