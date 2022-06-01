<script setup>
import { computed, reactive, ref, watch } from 'vue'
import { useSessionStore } from '../store/session'

import Chart from 'primevue/chart'

const session = useSessionStore()

const chart = ref(null)
const xAxisColumn = ref(null)
const yAxisColumn = ref(null)

const tooltip = reactive({
  x: null,
  y: null,
  items: []
})

const chartData = computed(() => {
  const data = session.frames.map(item => ({
    x: item[xAxisColumn.value],
    y: item[yAxisColumn.value]
  }))

  return { datasets: [{ data }] }
})

const tooltipHandler = ctx => {
  if (!ctx.tooltip.opacity) {
    tooltip.items = []
    return
  }

  const rect = ctx.chart.canvas.getBoundingClientRect()
  tooltip.x = rect.left + ctx.tooltip.caretX
  tooltip.y = rect.top + ctx.tooltip.caretY + 5
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
      radius: ctx => {
        return isSelected(ctx.dataIndex) ? 5 : 3
      }
    }
  },
  interaction: {
    mode: 'nearest',
    intersect: false
  },
  plugins: {
    legend: {
      display: false
    },
    title: {
      display: true,
      text: `${xAxisColumn.value} vs ${yAxisColumn.value}`
    },
    tooltip: {
      enabled: false,
      position: 'nearest',
      external: tooltipHandler
    }
  },
  scales: {
    x: {
      title: {
        display: true,
        text: xAxisColumn.value
      }
    },
    y: {
      title: {
        display: true,
        text: yAxisColumn.value
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
</script>

<template>
  <div class="mb-2">
    <div class="mx-2 inline-block">
      <div class="mb-2 text-sm text-left">X Axis</div>
      <Dropdown
        v-model="xAxisColumn"
        :options="session.numericColumns"
        class="w-10rem"
        placeholder="select column"
      />
    </div>

    <div class="mx-2 inline-block">
      <div class="mb-2 text-sm text-left">Y Axis</div>
      <Dropdown
        v-model="yAxisColumn"
        :options="session.numericColumns"
        class="w-10rem"
        placeholder="select column"
      />
    </div>
  </div>

  <Chart
    v-if="xAxisColumn && yAxisColumn"
    ref="chart"
    :data="chartData"
    :options="chartOptions"
    type="scatter"
    @click="onClick"
  />

  <Skeleton v-else class="mt-4" animation="none" height="14rem" />

  <div
    class="tooltip"
    :style="{
      top: `${tooltip.y}px`,
      left: `${tooltip.x}px`,
      opacity: tooltip.items.length ? 1 : 0
    }"
  >
    <div v-if="tooltip.items.length" class="p-2">
      <div>{{ xAxisColumn }}: {{ tooltip.items[0][xAxisColumn] }}</div>
      <div>{{ yAxisColumn }}: {{ tooltip.items[0][yAxisColumn] }}</div>
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
}
</style>
