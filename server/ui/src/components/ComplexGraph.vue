<script setup>
import { computed, reactive, ref, toRef, watch, watchEffect } from 'vue'
import { useSessionStore } from '../store/session'
import { useRadarData, useScatterData } from '../helpers/graphs'

import Chart from 'primevue/chart'

const props = defineProps({
  graph: Object,
  fullscreen: Boolean
})

const session = useSessionStore()
const graph = toRef(props, 'graph')

const chart = ref(null)
const settings = ref(null)

const radarData = useRadarData(session, graph)
const scatterData = useScatterData(session, graph)
const chartData = computed(() => {
  return graph.value.type === 'radar' ? radarData.value : scatterData.value
})

const tooltip = reactive({
  x: null,
  y: null,
  items: []
})

const fontMultiplier = computed(() => (props.fullscreen ? 1.2 : 1))
const radiusMultiplier = computed(() => (props.fullscreen ? 1.5 : 1))

const tooltipHandler = ctx => {
  if (!ctx.tooltip.opacity) {
    tooltip.items = []
    return
  }

  const rect = ctx.chart.canvas.getBoundingClientRect()
  tooltip.x = Math.min(rect.left + ctx.tooltip.caretX, window.innerWidth - 90)
  tooltip.y = rect.top + ctx.tooltip.caretY + 5 * radiusMultiplier.value
  tooltip.items = ctx.tooltip.dataPoints.map(d => {
    return session.frames.find(f => f.index === d.raw.index)
  })
}

const isSelected = ctx => {
  const frame = session.frames.find(f => f.index === ctx.raw.index)
  return session.selectedFrame === frame
}

const chartOptions = computed(() => ({
  elements: {
    point: {
      borderColor: '#fff',
      borderWidth: 0,
      // hoverBorderWidth: 1,
      backgroundColor: ctx => {
        return isSelected(ctx) ? '#fff' : '#fff4'
      },
      hoverRadius: () => 4 * radiusMultiplier.value,
      radius: ctx => {
        const r = isSelected(ctx) ? 4 : 3
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
      display: graph.value.type === 'scatter',
      font: { size: 16 * fontMultiplier.value, weight: 'normal' },
      text: `${graph.value.yColumn} vs ${graph.value.xColumn}`
    },
    tooltip: {
      enabled: false,
      position: 'nearest',
      filter: item => {
        return graph.value.type === 'scatter' ? item.datasetIndex === 0 : true
      },
      external: tooltipHandler
    }
  },
  scales: {
    r: {
      display: graph.value.type === 'radar',
      ticks: { display: false },
      suggestedMin: -0.2,
      suggestedMax: 1
    },
    x: {
      display: graph.value.type === 'scatter',
      ticks: { font: { size: 14 * fontMultiplier.value } },
      title: {
        display: true,
        font: { size: 14 * fontMultiplier.value },
        text: graph.value.xColumn
      }
    },
    y: {
      display: graph.value.type === 'scatter',
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

const roundValue = value => {
  return Math.round(value * 100) / 100
}
</script>

<template>
  <div class="mb-5">
    <Chart
      v-if="chartData"
      ref="chart"
      :key="graph.type"
      :data="chartData"
      :options="chartOptions"
      :type="graph.type"
      @click="onClick"
    />

    <Skeleton v-else height="auto" style="aspect-ratio: 2/1" animation="none" />
  </div>

  <div class="mb-2 flex justify-content-center gap-2">
    <template v-if="graph.type === 'scatter'">
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
    </template>

    <template v-else-if="graph.type === 'radar'">
      <span class="p-float-label">
        <MultiSelect
          v-model="graph.rColumns"
          :options="session.numericColumns"
          :max-selected-labels="0.1"
          placeholder="select columns"
          selected-items-label="{0} columns"
          class="w-8rem p-inputwrapper-filled"
        />
        <label>Axes</label>
      </span>

      <span class="p-float-label">
        <MultiSelect
          v-model="graph.frames"
          :max-selected-labels="-1"
          :options="session.frames"
          :option-label="f => `${f.index + 1} - ${f[session.nameColumn]}`"
          option-value="index"
          class="w-full p-inputwrapper-filled"
          placeholder="all frames"
          selected-items-label="{0} frames"
        />
        <label>Frames</label>
      </span>
    </template>

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
      <span class="mt-3 p-float-label">
        <Dropdown
          v-model="graph.type"
          :options="['radar', 'scatter']"
          class="w-full"
        />
        <label>Graph Type</label>
      </span>

      <span v-if="graph.type === 'scatter'" class="mt-3 p-float-label">
        <MultiSelect
          v-model="graph.frames"
          :max-selected-labels="-1"
          :options="session.frames"
          :option-label="f => `${f.index + 1} - ${f[session.nameColumn]}`"
          option-value="index"
          class="w-full p-inputwrapper-filled"
          placeholder="all frames"
          selected-items-label="{0} frames"
        />
        <label>Included Frames</label>
      </span>

      <template v-if="graph.type === 'scatter'">
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
      </template>
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
      <template v-if="graph.type === 'scatter'">
        <div class="flex justify-content-between">
          <div>{{ graph.xColumn }}:</div>
          <div>{{ roundValue(tooltip.items[0][graph.xColumn]) }}</div>
        </div>
        <div class="flex justify-content-between">
          <div>{{ graph.yColumn }}:</div>
          <div>{{ roundValue(tooltip.items[0][graph.yColumn]) }}</div>
        </div>
      </template>

      <template v-else-if="graph.type === 'radar'">
        <div v-for="col in graph.rColumns" class="flex justify-content-between">
          <div>{{ col }}:</div>
          <div>{{ roundValue(tooltip.items[0][col]) }}</div>
        </div>
      </template>
    </div>

    <div
      v-for="item in tooltip.items.slice(0, 3)"
      :key="item.index"
      :style="{
        background: session.selectedFrame === item ? '#fff2' : '#0000'
      }"
      class="p-2 text-center"
    >
      <img
        v-if="session.getImage(item.index)"
        :src="session.getImage(item.index)"
        class="h-4rem"
      />
      <div v-else class="h-4rem inline-flex align-items-center">
        <i class="pi pi-exclamation-triangle text-500 text-5xl" />
      </div>
      <div class="text-xs">
        {{ item.index + 1 }} - {{ item[session.nameColumn] }}
      </div>
    </div>

    <div v-if="tooltip.items.length > 3">
      ... {{ tooltip.items.length - 3 }} more
    </div>
  </div>
</template>

<style scoped>
.tooltip {
  position: absolute;
  width: 160px;
  color: #fff;
  background-color: #0008;
  border: 1px solid #4448;
  border-radius: 4px;
  pointer-events: none;
  transform: translate(-50%, 0);
  transition: all 0.3s ease;
  user-select: none;
  white-space: nowrap;
  z-index: 1000;
}
</style>
