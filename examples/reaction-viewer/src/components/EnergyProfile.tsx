import { useMemo } from 'react'

interface EnergyProfileProps {
  energies: number[]
  currentFrame: number
  width?: number
  height?: number
}

/** SVG-based energy vs frame mini-chart overlaid on the viewer */
export default function EnergyProfile({ energies, currentFrame, width = 280, height = 100 }: EnergyProfileProps) {
  const { pathD, maxE, minE, tsFrame, markerX, markerY } = useMemo(() => {
    const n = energies.length
    if (n === 0) return { pathD: '', maxE: 0, minE: 0, tsFrame: 0, markerX: 0, markerY: 0 }

    const pad = 8
    const w = width - pad * 2
    const h = height - pad * 2

    let mn = Infinity, mx = -Infinity, tsF = 0
    for (let i = 0; i < n; i++) {
      if (energies[i] < mn) mn = energies[i]
      if (energies[i] > mx) { mx = energies[i]; tsF = i }
    }
    const range = mx - mn || 1

    const points: string[] = []
    for (let i = 0; i < n; i++) {
      const x = pad + (i / (n - 1)) * w
      const y = pad + (1 - (energies[i] - mn) / range) * h
      points.push(`${i === 0 ? 'M' : 'L'}${x.toFixed(1)},${y.toFixed(1)}`)
    }

    const cx = pad + (currentFrame / (n - 1)) * w
    const cy = pad + (1 - (energies[currentFrame] - mn) / range) * h

    return { pathD: points.join(' '), maxE: mx, minE: mn, tsFrame: tsF, markerX: cx, markerY: cy }
  }, [energies, currentFrame, width, height])

  if (energies.length === 0) return null

  return (
    <div className="energy-profile">
      <svg width={width} height={height + 20} viewBox={`0 0 ${width} ${height + 20}`}>
        {/* Gradient fill under curve */}
        <defs>
          <linearGradient id="energyGrad" x1="0" y1="0" x2="0" y2="1">
            <stop offset="0%" stopColor="#ef4444" stopOpacity="0.4" />
            <stop offset="100%" stopColor="#3b82f6" stopOpacity="0.1" />
          </linearGradient>
        </defs>

        {/* Area fill */}
        <path
          d={`${pathD} L${width - 8},${height - 8} L8,${height - 8} Z`}
          fill="url(#energyGrad)"
        />

        {/* Energy curve */}
        <path d={pathD} fill="none" stroke="#f97316" strokeWidth="2" strokeLinejoin="round" />

        {/* Current frame marker */}
        <circle cx={markerX} cy={markerY} r="5" fill="#ef4444" stroke="#fff" strokeWidth="1.5" />

        {/* TS label */}
        <text x={width - 8} y={height + 14} fill="#888" fontSize="10" textAnchor="end">
          TS frame {tsFrame} • E_a = {(maxE - energies[0]).toFixed(1)}
        </text>
        <text x={8} y={height + 14} fill="#888" fontSize="10">
          Frame {currentFrame}/{energies.length - 1}
        </text>
      </svg>
    </div>
  )
}
