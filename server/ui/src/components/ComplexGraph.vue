<script setup>
import { computed, reactive, ref, watch } from 'vue'

import Chart from 'primevue/chart'

const props = defineProps({
  columns: Array,
  columnTypes: Object,
  complexIndex: Number,
  data: Array,
  images: Object,
  nameColumn: String,
  selectedFrame: Object
})

const emit = defineEmits(['update:selectedFrame'])

const chart = ref(null)
const xAxisColumn = ref(null)
const yAxisColumn = ref(null)

const tooltip = reactive({
  x: null,
  y: null,
  items: []
})

const numericColumns = computed(() => {
  return props.columns.filter(c => props.columnTypes[c] === 'numeric')
})

const data = computed(() => {
  const data = props.data.map(item => ({
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
  tooltip.items = ctx.tooltip.dataPoints.map(d => props.data[d.dataIndex])
}

const isSelected = index => {
  return props.selectedFrame === props.data[index]
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
    title: {
      display: true,
      text: `${xAxisColumn.value} vs ${yAxisColumn.value}`
    },
    legend: {
      display: false
    },
    tooltip: {
      enabled: false,
      position: 'nearest',
      external: tooltipHandler
    }
  }
}))

watch(
  () => props.selectedFrame,
  () => chart.value?.refresh(),
  { deep: true }
)

const onClick = () => {
  if (!tooltip.items.length) return
  let index = tooltip.items.indexOf(props.selectedFrame)
  index = index === -1 ? 0 : (index + 1) % tooltip.items.length
  emit('update:selectedFrame', tooltip.items[index])
}
</script>

<template>
  <div class="mx-2 inline-block">
    <div class="mb-2 text-sm text-left">X Axis</div>
    <Dropdown
      v-model="xAxisColumn"
      :options="numericColumns"
      class="w-15rem"
      placeholder="select column"
    />
  </div>

  <div class="mx-2 inline-block">
    <div class="mb-2 text-sm text-left">Y Axis</div>
    <Dropdown
      v-model="yAxisColumn"
      :options="numericColumns"
      class="w-15rem"
      placeholder="select column"
    />
  </div>

  <Chart
    v-if="xAxisColumn && yAxisColumn"
    ref="chart"
    :data="data"
    :options="chartOptions"
    type="scatter"
    @click="onClick"
  />

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
      :style="{ background: selectedFrame === item ? '#fff2' : '#0000' }"
      class="p-2 text-center"
    >
      <img
        :src="
          'data:image/png;base64,' +
          props.images[props.complexIndex + '-' + item.index]
        "
        class="h-4rem"
      />
      <div>{{ item[props.nameColumn] }}</div>
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
