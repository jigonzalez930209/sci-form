import { Canvas } from '@react-three/fiber'
import { OrbitControls, Environment, ContactShadows } from '@react-three/drei'
import type { MoleculeFrame } from '../types'
import Molecule from './Molecule'

interface MoleculeViewerProps {
  frame: MoleculeFrame | null
}

/** Full 3D viewer canvas with lighting, controls, and molecule rendering */
export default function MoleculeViewer({ frame }: MoleculeViewerProps) {
  return (
    <div className="viewer-canvas">
      <Canvas
        camera={{ position: [0, 0, 8], fov: 50, near: 0.1, far: 100 }}
        gl={{ antialias: true, alpha: true }}
        dpr={[1, 2]}
      >
        {/* Lighting */}
        <ambientLight intensity={0.5} />
        <directionalLight position={[5, 5, 5]} intensity={1.0} castShadow />
        <directionalLight position={[-3, 3, -3]} intensity={0.4} />
        <pointLight position={[0, -3, 2]} intensity={0.3} color="#4d9fff" />

        {/* Environment for reflections */}
        <Environment preset="studio" />

        {/* Subtle contact shadows */}
        <ContactShadows
          position={[0, -3, 0]}
          opacity={0.3}
          scale={10}
          blur={2}
        />

        {/* Molecule */}
        {frame && <Molecule frame={frame} />}

        {/* Orbit controls */}
        <OrbitControls
          enableDamping
          dampingFactor={0.08}
          minDistance={2}
          maxDistance={30}
          enablePan
        />
      </Canvas>
    </div>
  )
}
