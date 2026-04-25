import { useState, useRef, useCallback, useEffect } from 'react'

interface UseAnimationOptions {
  totalFrames: number
  fps?: number
  speed?: number
  loop?: boolean
}

export function useAnimation({ totalFrames, fps = 12, speed = 1, loop = true }: UseAnimationOptions) {
  const [currentFrame, setCurrentFrame] = useState(0)
  const [isPlaying, setIsPlaying] = useState(false)
  const rafRef = useRef<number>(0)
  const lastTimeRef = useRef(0)

  const play = useCallback(() => setIsPlaying(true), [])
  const pause = useCallback(() => setIsPlaying(false), [])
  const toggle = useCallback(() => setIsPlaying((p) => !p), [])
  const goTo = useCallback((f: number) => {
    setCurrentFrame(Math.max(0, Math.min(totalFrames - 1, f)))
  }, [totalFrames])
  const reset = useCallback(() => { setCurrentFrame(0); setIsPlaying(false) }, [])

  useEffect(() => {
    if (!isPlaying || totalFrames <= 1) return

    const interval = 1000 / (fps * speed)
    lastTimeRef.current = performance.now()

    const tick = (now: number) => {
      const elapsed = now - lastTimeRef.current
      if (elapsed >= interval) {
        lastTimeRef.current = now - (elapsed % interval)
        setCurrentFrame((prev) => {
          const next = prev + 1
          if (next >= totalFrames) {
            if (loop) return 0
            setIsPlaying(false)
            return prev
          }
          return next
        })
      }
      rafRef.current = requestAnimationFrame(tick)
    }

    rafRef.current = requestAnimationFrame(tick)
    return () => cancelAnimationFrame(rafRef.current)
  }, [isPlaying, totalFrames, fps, speed, loop])

  // Reset frame when totalFrames changes
  useEffect(() => {
    setCurrentFrame((prev) => Math.min(prev, totalFrames - 1))
  }, [totalFrames])

  return { currentFrame, isPlaying, play, pause, toggle, goTo, reset }
}
