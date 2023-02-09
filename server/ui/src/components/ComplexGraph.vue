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
const isRadar = computed(() => graph.value.type === 'radar')
const isScatter = computed(() => graph.value.type === 'scatter')

const chart = ref(null)
const settings = ref(null)

const radarData = useRadarData(session, graph)
const scatterData = useScatterData(session, graph)
const chartData = computed(() => {
  return isRadar.value ? radarData.value : scatterData.value
})

const tooltip = reactive({
  label: null,
  x: null,
  y: null,
  items: [],
  start: null,
  total: null,
  nextSelection: null
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

  const total = ctx.tooltip.dataPoints.length
  const selectedIndex = ctx.tooltip.dataPoints.findIndex(d => {
    return d.raw.id === session.selectedFrame?.id
  })

  if (total > 0) {
    if (isRadar.value) {
      tooltip.label = ctx.tooltip.dataPoints[0].label
    }

    const nextIndex = (selectedIndex + 1) % total
    tooltip.nextSelection = ctx.tooltip.dataPoints[nextIndex].raw.id
  }

  const start = Math.floor((selectedIndex === -1 ? 0 : selectedIndex) / 3)
  let items = ctx.tooltip.dataPoints.slice(start * 3, start * 3 + 3).map(d => {
    return session.frames.find(f => f.id === d.raw.id)
  })

  tooltip.items = items
  tooltip.start = start
  tooltip.total = total
}

const isSelected = ctx => {
  const id = ctx.raw ? ctx.raw.id : ctx.dataset.data[0].id
  const frame = session.frames.find(f => f.id === id)
  return session.selectedFrame === frame
}

const chartOptions = computed(() => ({
  elements: {
    line: {
      borderColor: ctx => {
        return isSelected(ctx) ? '#fff' : '#fff2'
      },
      borderJoinStyle: 'round',
      borderWidth: () => 4 * radiusMultiplier.value,
      fill: false
    },
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
      display: isScatter.value,
      font: { size: 16 * fontMultiplier.value, weight: 'normal' },
      text: `${graph.value.yColumn} vs ${graph.value.xColumn}`
    },
    tooltip: {
      enabled: false,
      position: 'nearest',
      filter: item => {
        return isScatter.value ? item.datasetIndex === 0 : true
      },
      external: tooltipHandler
    }
  },
  scales: {
    r: {
      display: isRadar.value,
      min: -0.2,
      max: 1,
      pointLabels: {
        font: { size: 12 * fontMultiplier.value }
      },
      ticks: { display: false }
    },
    x: {
      display: isScatter.value,
      ticks: { font: { size: 14 * fontMultiplier.value } },
      title: {
        display: true,
        font: { size: 14 * fontMultiplier.value },
        text: graph.value.xColumn
      }
    },
    y: {
      display: isScatter.value,
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

const cycleSelection = () => {
  if (!tooltip.items.length) return
  session.selectFrame(tooltip.nextSelection)
}

const syncSelectionWithTable = to => {
  if (to) {
    const ids = graph.value.frames
    session.selectedFrames = ids.map(id => session.find(f => f.id === id))
  } else {
    const ids = session.selectedFrames.map(f => f.id)
    graph.value.frames = ids
  }
  session.selectionMode = true
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
  <div
    :class="fullscreen ? 'h-full flex-center' : 'flex-column'"
    class="flex gap-5"
  >
    <div
      class="flex-grow-1"
      :style="graph.type === 'radar' ? 'max-width: calc(100vh - 80px)' : ''"
    >
      <Chart
        v-if="chartData"
        ref="chart"
        :key="graph.type"
        :data="chartData"
        :options="chartOptions"
        :type="graph.type"
        @click="cycleSelection"
      />

      <Skeleton
        v-else
        :style="{ 'aspect-ratio': isScatter ? '2/1' : '1/1' }"
        animation="none"
        height="auto"
      />
    </div>

    <!-- controls -->
    <div
      :class="{ 'w-20rem flex-column mt-4': fullscreen }"
      class="mb-2 flex gap-1"
    >
      <template v-if="isScatter">
        <span
          :class="fullscreen ? 'mt-4' : 'flex-grow-1'"
          class="p-float-label"
        >
          <Dropdown
            v-model="graph.yColumn"
            :options="session.numericColumns"
            class="w-full"
          />
          <label>Y Axis</label>
        </span>

        <div
          v-if="fullscreen || graph.xColumn || graph.yColumn"
          class="text-center"
        >
          <Button
            v-tooltip.bottom="'swap axes'"
            :class="{ 'rotate-90': !fullscreen }"
            class="p-button-secondary p-button-text"
            icon="pi pi-sort-alt"
            @click="swapAxes"
          />
        </div>

        <span :class="fullscreen ? '' : 'flex-grow-1'" class="p-float-label">
          <Dropdown
            v-model="graph.xColumn"
            :options="session.numericColumns"
            class="w-full"
          />
          <label>X Axis</label>
        </span>
      </template>

      <template v-else-if="isRadar">
        <span
          :class="fullscreen ? 'mt-4' : 'flex-grow-1'"
          class="p-float-label"
        >
          <MultiSelect
            v-model="graph.rColumns"
            :options="session.numericColumns"
            :max-selected-labels="0.1"
            class="w-full p-inputwrapper-filled"
            placeholder="select columns"
            scroll-height="500px"
            selected-items-label="{0} columns"
          />
          <label>Axes</label>
        </span>

        <span
          :class="fullscreen ? 'mt-4' : 'flex-grow-1'"
          class="p-float-label"
        >
          <MultiSelect
            v-model="graph.frames"
            :max-selected-labels="-1"
            :options="session.frames"
            option-label="frame"
            option-value="id"
            class="w-full p-inputwrapper-filled"
            placeholder="all frames"
            scroll-height="500px"
            selected-items-label="{0} frames"
          />
          <label>Frames</label>
        </span>
      </template>

      <div class="flex flex-grow-0">
        <Button
          v-tooltip.bottom="'settings'"
          :label="fullscreen ? 'settings' : ''"
          class="mx-auto p-button-secondary p-button-text"
          icon="pi pi-cog"
          @click="e => settings.toggle(e)"
        />

        <template v-if="!fullscreen">
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
    </div>
  </div>

  <Divider v-if="!fullscreen" />

  <!-- settings -->
  <OverlayPanel ref="settings">
    <div class="mt-3 w-12rem flex flex-column gap-5">
      <span class="p-float-label">
        <Dropdown
          v-model="graph.type"
          :options="['radar', 'scatter']"
          class="w-full"
          @change="settings.hide()"
        />
        <label>Graph Type</label>
      </span>

      <template v-if="isScatter">
        <span class="p-float-label">
          <MultiSelect
            v-model="graph.frames"
            :max-selected-labels="-1"
            :options="session.frames"
            option-label="frame"
            option-value="id"
            class="w-full p-inputwrapper-filled"
            placeholder="all frames"
            scroll-height="500px"
            selected-items-label="{0} frames"
          />
          <label>Frames</label>
        </span>

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

      <span class="w-full p-float-label">
        <div class="flex p-buttonset p-inputwrapper-filled">
          <Button
            :label="`${graph.frames.length}`"
            class="w-full p-button-secondary p-button-outlined"
            icon="pi pi-arrow-left"
            @click="syncSelectionWithTable(true)"
          />
          <Button
            :label="`${session.selectedFrames.length}`"
            class="w-full p-button-secondary p-button-outlined"
            icon="pi pi-arrow-right"
            @click="syncSelectionWithTable(false)"
          />
        </div>
        <label>Sync selection with table</label>
      </span>
    </div>
  </OverlayPanel>

  <!-- tooltip -->
  <div
    class="tooltip"
    :style="{
      top: `${tooltip.y}px`,
      left: `${tooltip.x}px`,
      opacity: tooltip.items.length ? 1 : 0
    }"
  >
    <div v-if="tooltip.items.length" class="p-2">
      <template v-if="isScatter">
        <div class="flex justify-content-between">
          <div>{{ graph.xColumn }}:</div>
          <div>{{ roundValue(tooltip.items[0][graph.xColumn]) }}</div>
        </div>
        <div class="flex justify-content-between">
          <div>{{ graph.yColumn }}:</div>
          <div>{{ roundValue(tooltip.items[0][graph.yColumn]) }}</div>
        </div>
      </template>

      <template v-else-if="isRadar">
        <div
          v-for="col in graph.rColumns"
          :class="{ 'font-bold': col === tooltip.label }"
          class="flex justify-content-between"
        >
          <div>{{ col }}:</div>
          <div>{{ roundValue(tooltip.items[0][col]) }}</div>
        </div>
      </template>
    </div>

    <div
      v-for="item in tooltip.items"
      :key="item.id"
      :style="{
        background: session.selectedFrame === item ? '#fff2' : '#0000'
      }"
      class="p-2 text-center"
    >
      <img
        v-if="session.getImage(item.id)"
        :src="session.getImage(item.id)"
        class="h-4rem"
      />
      <div v-else class="h-4rem inline-flex align-items-center">
        <i class="pi pi-exclamation-triangle text-500 text-5xl" />
      </div>
      <div class="text-xs">
        {{ item.frame }}
      </div>
    </div>

    <div v-if="tooltip.total > 3">{{ tooltip.total }} frames</div>
  </div>
</template>

<style scoped>
.tooltip {
  position: absolute;
  width: 160px;
  color: #fff;
  background-color: #0008;
  backdrop-filter: blur(4px);
  border-radius: 4px;
  pointer-events: none;
  transform: translate(-50%, 0);
  transition: all 0.3s ease;
  user-select: none;
  white-space: nowrap;
  z-index: 1000;
}
</style>
