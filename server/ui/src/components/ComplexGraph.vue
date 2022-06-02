<script setup>
import { computed, reactive, ref, toRef, watch } from 'vue'
import { useSessionStore } from '../store/session'

import Chart from 'primevue/chart'

const props = defineProps({
  graph: Object,
  fullscreen: Boolean
})

const session = useSessionStore()
const graph = toRef(props, 'graph')

const chart = ref(null)
const tooltip = reactive({
  x: null,
  y: null,
  items: []
})

const chartData = computed(() => {
  const data = session.frames.map(item => ({
    x: item[graph.value.xColumn],
    y: item[graph.value.yColumn]
  }))

  return { datasets: [{ data }] }
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
      type="scatter"
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
