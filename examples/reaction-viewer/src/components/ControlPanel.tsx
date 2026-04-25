import type { Algorithm, ReactionDefinition } from '../types'

interface ControlPanelProps {
  reactions: ReactionDefinition[]
  selectedReaction: string
  onSelectReaction: (id: string) => void
  algorithm: Algorithm
  onAlgorithmChange: (a: Algorithm) => void
  nFrames: number
  onNFramesChange: (n: number) => void
  currentFrame: number
  onFrameChange: (f: number) => void
  totalFrames: number
  isPlaying: boolean
  onPlayPause: () => void
  speed: number
  onSpeedChange: (s: number) => void
  onRegenerate: () => void
}

const ALGORITHMS: { value: Algorithm; label: string }[] = [
  { value: 'uff', label: 'UFF' },
  { value: 'pm3', label: 'PM3' },
  { value: 'xtb', label: 'GFN0-xTB' },
  { value: 'gfn1', label: 'GFN1-xTB' },
  { value: 'gfn2', label: 'GFN2-xTB' },
  { value: 'hf3c', label: 'HF-3c' },
]

const FRAME_OPTIONS = [11, 15, 21, 31, 41, 51]

export default function ControlPanel({
  reactions,
  selectedReaction,
  onSelectReaction,
  algorithm,
  onAlgorithmChange,
  nFrames,
  onNFramesChange,
  currentFrame,
  onFrameChange,
  totalFrames,
  isPlaying,
  onPlayPause,
  speed,
  onSpeedChange,
  onRegenerate,
}: ControlPanelProps) {
  return (
    <div className="control-panel">
      {/* ── Reaction selector ────────────── */}
      <div className="control-group">
        <label>Reaction</label>
        <div className="reaction-tabs">
          {reactions.map((r) => (
            <button
              key={r.id}
              className={`tab ${selectedReaction === r.id ? 'active' : ''}`}
              onClick={() => onSelectReaction(r.id)}
              title={r.description}
            >
              {r.name}
            </button>
          ))}
        </div>
      </div>

      {/* ── Algorithm + Frames ───────────── */}
      <div className="control-row">
        <div className="control-group">
          <label>Algorithm</label>
          <select
            value={algorithm}
            onChange={(e) => onAlgorithmChange(e.target.value as Algorithm)}
          >
            {ALGORITHMS.map((a) => (
              <option key={a.value} value={a.value}>{a.label}</option>
            ))}
          </select>
        </div>

        <div className="control-group">
          <label>Frames</label>
          <select
            value={nFrames}
            onChange={(e) => onNFramesChange(Number(e.target.value))}
          >
            {FRAME_OPTIONS.map((n) => (
              <option key={n} value={n}>{n}</option>
            ))}
          </select>
        </div>

        <div className="control-group">
          <label>Speed</label>
          <select value={speed} onChange={(e) => onSpeedChange(Number(e.target.value))}>
            <option value={0.25}>0.25×</option>
            <option value={0.5}>0.5×</option>
            <option value={1}>1×</option>
            <option value={2}>2×</option>
            <option value={4}>4×</option>
          </select>
        </div>

        <button className="btn btn-accent" onClick={onRegenerate} title="Regenerate frames with selected algorithm">
          ⟳ Generate
        </button>
      </div>

      {/* ── Playback controls ────────────── */}
      <div className="control-row playback">
        <button className="btn btn-play" onClick={onPlayPause}>
          {isPlaying ? '⏸' : '▶'}
        </button>

        <input
          type="range"
          min={0}
          max={totalFrames - 1}
          value={currentFrame}
          onChange={(e) => onFrameChange(Number(e.target.value))}
          className="frame-slider"
        />

        <span className="frame-label">
          {currentFrame} / {totalFrames - 1}
        </span>
      </div>
    </div>
  )
}
