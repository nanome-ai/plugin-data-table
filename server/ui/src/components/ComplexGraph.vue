<script setup>
import { computed, reactive, ref, toRef, watch, watchEffect } from 'vue'
import regression from 'regression'
import { useSessionStore } from '../store/session'

import Chart from 'primevue/chart'

const props = defineProps({
  graph: Object,
  fullscreen: Boolean
})

const session = useSessionStore()
const graph = toRef(props, 'graph')

const chart = ref(null)
const settings = ref(null)
const chartData = ref([])

const tooltip = reactive({
  x: null,
  y: null,
  items: []
})

watchEffect(() => {
  const data = session.frames.map(item => ({
    x: +item[graph.value.xColumn],
    y: +item[graph.value.yColumn]
  }))

  const datasets = [{ data }]
  chartData.value = { datasets }

  graph.value.reg.error = false

  if (graph.value.type === 'scatter' && graph.value.reg.type !== 'none') {
    const points = data
      .map(point => [point.x, point.y])
      .sort((a, b) => a[0] - b[0])
    const r = regression[graph.value.reg.type](points, {
      order: graph.value.reg.order
    })

    graph.value.reg.error = isNaN(r.r2)
    graph.value.reg.eq = r.string.replace(/\^([-.\d()x]+)/g, '<sup>$1</sup>')
    graph.value.reg.r2 = r.r2
    if (graph.value.reg.error) return

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
      borderColor: '#' + graph.value.reg.color,
      pointRadius: 0,
      pointHitRadius: 0,
      pointHoverRadius: 0
    })
  }
})

const fontMultiplier = computed(() => (props.fullscreen ? 1.2 : 1))
const radiusMultiplier = computed(() => (props.fullscreen ? 1.5 : 1))

const tooltipHandler = ctx => {
  if (!ctx.tooltip.opacity) {
    tooltip.items = []
    return
  }

  const rect = ctx.chart.canvas.getBoundingClientRect()
  tooltip.x = rect.left + ctx.tooltip.caretX
  tooltip.y = rect.top + ctx.tooltip.caretY + 5 * radiusMultiplier.value
  tooltip.items = ctx.tooltip.dataPoints.map(d => session.frames[d.dataIndex])
}

const isSelected = index => {
  return session.selectedFrame === session.frames[index]
}

const chartOptions = computed(() => ({
  elements: {
    point: {
      borderColor: '#fff',
      borderWidth: 0,
      // hoverBorderWidth: 1,
      backgroundColor: ctx => {
        return isSelected(ctx.dataIndex) ? '#fff' : '#fff4'
      },
      hoverRadius: () => 4 * radiusMultiplier.value,
      radius: ctx => {
        const r = isSelected(ctx.dataIndex) ? 4 : 3
        return r * radiusMultiplier.value
      }
    }
  },
  interaction: {
    mode: 'nearest',
    intersect: false
  },
  plugins: {
    legend: { display: false },
    title: {
      display: true,
      font: { size: 16 * fontMultiplier.value, weight: 'normal' },
      text: `${graph.value.yColumn} vs ${graph.value.xColumn}`
    },
    tooltip: {
      enabled: false,
      position: 'nearest',
      filter: item => item.datasetIndex === 0,
      external: tooltipHandler
    }
  },
  scales: {
    x: {
      ticks: { font: { size: 14 * fontMultiplier.value } },
      title: {
        display: true,
        font: { size: 14 * fontMultiplier.value },
        text: graph.value.xColumn
      }
    },
    y: {
      ticks: { font: { size: 14 * fontMultiplier.value } },
      title: {
        display: true,
        font: { size: 14 * fontMultiplier.value },
        text: graph.value.yColumn
      }
    }
  }
}))

watch(
  () => session.selectedFrame,
  () => chart.value?.refresh(),
  { deep: true }
)

const onClick = () => {
  if (!tooltip.items.length) return
  let index = tooltip.items.indexOf(session.selectedFrame)
  index = index === -1 ? 0 : (index + 1) % tooltip.items.length
  session.selectFrame(tooltip.items[index].index)
}

const swapAxes = () => {
  const x = graph.value.xColumn
  graph.value.xColumn = graph.value.yColumn
  graph.value.yColumn = x
}
</script>

<template>
  <div class="mb-5">
    <Chart
      v-if="graph.xColumn && graph.yColumn"
      ref="chart"
      :data="chartData"
      :options="chartOptions"
      :type="graph.type"
      @click="onClick"
    />

    <Skeleton v-else height="auto" style="aspect-ratio: 2/1" animation="none" />
  </div>

  <div class="mb-2 flex justify-content-center gap-2">
    <span class="p-float-label">
      <Dropdown
        v-model="graph.yColumn"
        :options="session.numericColumns"
        class="w-8rem"
      />
      <label>Y Axis</label>
    </span>

    <Button
      v-tooltip.bottom="'swap axes'"
      class="p-button-secondary p-button-text rotate-90"
      icon="pi pi-sort-alt"
      @click="swapAxes"
    />

    <span class="p-float-label">
      <Dropdown
        v-model="graph.xColumn"
        :options="session.numericColumns"
        class="w-8rem"
      />
      <label>X Axis</label>
    </span>

    <Button
      v-tooltip.bottom="'settings'"
      class="p-button-secondary p-button-text"
      icon="pi pi-cog"
      @click="e => settings.toggle(e)"
    />

    <template v-if="!props.fullscreen">
      <Button
        v-tooltip.bottom="'fullscreen'"
        class="p-button-secondary p-button-text"
        icon="pi pi-window-maximize"
        @click="session.selectGraph(graph)"
      />

      <Button
        v-tooltip.bottom="'remove graph'"
        class="p-button-danger p-button-text"
        icon="pi pi-trash"
        @click="session.removeGraph(graph)"
      />
    </template>
  </div>

  <Divider v-if="!props.fullscreen" />

  <OverlayPanel ref="settings">
    <div class="w-12rem flex flex-column gap-5">
      <!-- <span class="mt-3 p-float-label">
        <Dropdown
          v-model="graph.type"
          :options="['radar', 'scatter']"
          class="w-full"
        />
        <label>Graph Type</label>
      </span> -->

      <div v-if="graph.reg.type !== 'none'" class="text-center text-xs">
        <template v-if="!graph.reg.error">
          <div class="mb-2">r<sup>2</sup> = {{ graph.reg.r2 }}</div>
          <div v-html="graph.reg.eq" />
        </template>

        <span v-else class="p-error">error calculating regression</span>
      </div>

      <div class="mt-3 flex-center gap-2">
        <span class="p-float-label flex-grow-1">
          <Dropdown
            v-model="graph.reg.type"
            :class="{ 'p-invalid': graph.reg.error }"
            :options="[
              'none',
              'exponential',
              'linear',
              'logarithmic',
              'polynomial',
              'power'
            ]"
            class="w-full"
          />
          <label>Regression</label>
        </span>

        <ColorPicker v-model="graph.reg.color" />
      </div>

      <span v-if="graph.reg.type === 'polynomial'" class="p-float-label">
        <Dropdown
          v-model="graph.reg.order"
          :options="[2, 3, 4, 5, 6]"
          class="w-full"
        />
        <label>Polynomial Order</label>
      </span>
    </div>
  </OverlayPanel>

  <div
    class="tooltip"
    :style="{
      top: `${tooltip.y}px`,
      left: `${tooltip.x}px`,
      opacity: tooltip.items.length ? 1 : 0
    }"
  >
    <div v-if="tooltip.items.length" class="p-2">
      <div>{{ graph.xColumn }}: {{ tooltip.items[0][graph.xColumn] }}</div>
      <div>{{ graph.yColumn }}: {{ tooltip.items[0][graph.yColumn] }}</div>
    </div>
    <div
      v-for="item in tooltip.items.slice(0, 3)"
      :key="item.index"
      :style="{
        background: session.selectedFrame === item ? '#fff2' : '#0000'
      }"
      class="p-2 text-center"
    >
      <img :src="session.getImage(item.index)" class="h-4rem" />
      <div>{{ item[session.nameColumn] }}</div>
    </div>
    <div v-if="tooltip.items.length > 3">
      ... {{ tooltip.items.length - 3 }} more
    </div>
  </div>
</template>

<style scoped>
.tooltip {
  position: absolute;
  color: #fff;
  background-color: #0008;
  border: 1px solid #4448;
  border-radius: 4px;
  pointer-events: none;
  transform: translate(-50%, 0);
  transition: all 0.3s ease;
  user-select: none;
  z-index: 1000;
}
</style>
