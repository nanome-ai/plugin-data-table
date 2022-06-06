import { ref, watchEffect } from 'vue'
import regression from 'regression'

export const useRadarData = (session, graph) => {
  const graphData = ref(null)

  watchEffect(() => {
    const g = graph.value
    const cols = g.rColumns
    if (g.type !== 'radar' || !cols.length) {
      graphData.value = null
      return
    }

    const frames = session.frames.filter(
      f => !g.frames.length || g.frames.includes(f.index)
    )
    const datasets = []
    graphData.value = { labels: cols, datasets }

    const minValues = cols.map(c => Math.min(...frames.map(f => f[c])))
    const maxValues = cols.map(c => Math.max(...frames.map(f => f[c])))
    const diffValues = cols.map((c, i) => maxValues[i] - minValues[i])

    for (const frame of frames) {
      // remap values 0 to 1
      const data = cols.map((c, i) => {
        return (frame[c] - minValues[i]) / diffValues[i]
      })
      datasets.push({ data })
    }
  })

  return graphData
}

export const useScatterData = (session, graph) => {
  const graphData = ref(null)

  watchEffect(() => {
    const g = graph.value
    if (g.type !== 'scatter' || !g.xColumn || !g.yColumn) {
      graphData.value = null
      return
    }

    const data = session.frames
      .filter(f => !g.frames.length || g.frames.includes(f.index))
      .map(item => ({ x: +item[g.xColumn], y: +item[g.yColumn] }))

    const datasets = [{ data }]
    graphData.value = { datasets }

    g.reg.error = false

    if (g.reg.type !== 'none') {
      const points = data
        .map(point => [point.x, point.y])
        .sort((a, b) => a[0] - b[0])
      const r = regression[g.reg.type](points, {
        order: g.reg.order
      })

      g.reg.error = isNaN(r.r2)
      g.reg.eq = r.string.replace(/\^([-.\d()x]+)/g, '<sup>$1</sup>')
      g.reg.r2 = r.r2
      if (g.reg.error) return

      const xValues = points.map(point => point[0])
      const minX = Math.min(...xValues)
      const maxX = Math.max(...xValues)
      const step = (maxX - minX) / 300

      const rData = Array.from({ length: 300 }, (_, i) => ({
        x: minX + i * step,
        y: r.predict(minX + i * step)[1]
      }))

      const yValues = points.map(point => point[1])
      const minY = Math.min(...yValues)
      const maxY = Math.max(...yValues)
      const yPad = (maxY - minY) / 10

      // hide regression that goes out of bounds
      for (const p of rData) {
        if (p.y < minY - yPad || p.y > maxY + yPad) p.y = null
      }

      datasets.push({
        type: 'line',
        data: rData,
        borderColor: '#' + g.reg.color,
        pointRadius: 0,
        pointHitRadius: 0,
        pointHoverRadius: 0
      })
    }
  })

  return graphData
}
