import { useState, useMemo, useCallback } from 'react'
import MoleculeViewer from './components/MoleculeViewer'
import ControlPanel from './components/ControlPanel'
import EnergyProfile from './components/EnergyProfile'
import { useAnimation } from './hooks/useAnimation'
import { getReactions } from './data/reactions'
import type { Algorithm, ReactionDefinition } from './types'

export default function App() {
  const [nFrames, setNFrames] = useState(21)
  const [algorithm, setAlgorithm] = useState<Algorithm>('uff')
  const [selectedId, setSelectedId] = useState('h2-dissociation')
  const [speed, setSpeed] = useState(1)

  // Generate reactions with current frame count
  const reactions = useMemo(() => getReactions(nFrames), [nFrames])
  const reaction = useMemo(
    () => reactions.find((r) => r.id === selectedId) ?? reactions[0],
    [reactions, selectedId],
  )

  const { currentFrame, isPlaying, toggle, goTo } = useAnimation({
    totalFrames: reaction.frames.length,
    speed,
  })

  const frame = reaction.frames[currentFrame] ?? null
  const energies = useMemo(() => reaction.frames.map((f) => f.energy), [reaction])

  const handleSelectReaction = useCallback((id: string) => {
    setSelectedId(id)
    goTo(0)
  }, [goTo])

  const handleRegenerate = useCallback(() => {
    // In a real integration, this would call WASM NEB
    // For now, regenerate with new frame count (already reactive via useMemo)
    goTo(0)
  }, [goTo])

  return (
    <div className="app">
      <header className="app-header">
        <h1>sci-form — 3D Reaction Viewer</h1>
        <p className="subtitle">Interactive visualization of chemical reaction pathways</p>
      </header>

      <div className="viewer-container">
        {/* 3D Canvas */}
        <MoleculeViewer frame={frame} />

        {/* Energy profile overlay */}
        <EnergyProfile
          energies={energies}
          currentFrame={currentFrame}
        />

        {/* Reaction info overlay */}
        <div className="reaction-info">
          <h3>{reaction.name}</h3>
          <p>{reaction.description}</p>
          <code className="smirks">{reaction.smirks}</code>
        </div>
      </div>

      {/* Controls below viewer */}
      <ControlPanel
        reactions={reactions}
        selectedReaction={selectedId}
        onSelectReaction={handleSelectReaction}
        algorithm={algorithm}
        onAlgorithmChange={setAlgorithm}
        nFrames={nFrames}
        onNFramesChange={setNFrames}
        currentFrame={currentFrame}
        onFrameChange={goTo}
        totalFrames={reaction.frames.length}
        isPlaying={isPlaying}
        onPlayPause={toggle}
        speed={speed}
        onSpeedChange={setSpeed}
        onRegenerate={handleRegenerate}
      />
    </div>
  )
}
