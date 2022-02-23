<script setup>
import { computed, ref, watch } from 'vue'
import { useMux, useVModel } from '../composables.js'

import { FilterMatchMode, FilterOperator } from 'primevue/api'
import DataTable from 'primevue/datatable'
import Column from 'primevue/column'
import Image from 'primevue/image'

const OPERATOR_MAP = {
  [FilterMatchMode.EQUALS]: '=',
  [FilterMatchMode.NOT_EQUALS]: '≠',
  [FilterMatchMode.LESS_THAN]: '<',
  [FilterMatchMode.LESS_THAN_OR_EQUAL_TO]: '≤',
  [FilterMatchMode.GREATER_THAN]: '>',
  [FilterMatchMode.GREATER_THAN_OR_EQUAL_TO]: '≥'
}

const props = defineProps({
  columns: Array,
  columnTypes: Object,
  complexIndex: Number,
  hiddenFrames: Array,
  images: Object,
  loading: Boolean,
  modelValue: Array,
  multiSelect: Boolean,
  nameColumn: String,
  reorderable: Boolean,
  selectedFrame: Object,
  selectedRows: Array
})

const emit = defineEmits([
  'update:modelValue',
  'update:selectedFrame',
  'update:selectedRows'
])

const data = useVModel(props, emit)
const selectedFrame = useVModel(props, emit, 'selectedFrame')
const selectedRows = useVModel(props, emit, 'selectedRows')

const sortField = ref(null)
const selection = useMux(() => props.multiSelect, selectedRows, selectedFrame)

const filteredColumns = computed(() => {
  const columns = props.columns.filter(c => c !== props.nameColumn)
  columns.unshift(props.nameColumn)
  return columns
})

const filteredRows = computed(() => {
  return data.value.filter(row => !props.hiddenFrames.includes(row.index))
})

const filters = ref({})
const activeFilters = ref([])

watch(
  filters,
  () => {
    activeFilters.value = []
    filteredColumns.value.forEach(column => {
      const { constraints } = filters.value[column] || {}
      if (!constraints) return

      Object.entries(constraints).forEach(([index, constraint]) => {
        if (constraint.value === null) return
        let op =
          OPERATOR_MAP[constraint.matchMode] ||
          constraint.matchMode.replace(/[A-Z]/g, ' $&').trim().toLowerCase()
        const text = `${column} ${op} ${constraint.value}`
        activeFilters.value.push({ column, index, constraint, text })
      })
    })
  },
  { deep: true }
)

const resetFilters = () => {
  filters.value = {}
  for (const column of props.columns) {
    const matchMode = {
      text: FilterMatchMode.CONTAINS,
      numeric: FilterMatchMode.LESS_THAN
    }[props.columnTypes[column]]

    filters.value[column] = {
      operator: FilterOperator.AND,
      constraints: [{ matchMode, value: null }]
    }
  }
}

watch(() => props.columnTypes, resetFilters, { immediate: true })

const getImage = index => {
  const id = props.complexIndex + '-' + index
  return props.images[id]
}

const onColumnReorder = (...args) => {
  console.log('onColumnReorder', ...args)
}

const onRowReorder = e => {
  sortField.value = null
  data.value = e.value
}

const resetFilter = filter => {
  const constraints = filters.value[filter.column].constraints
  if (constraints.length > 1) {
    constraints.splice(filter.index, 1)
  } else {
    constraints[filter.index].value = null
  }
}
</script>

<template>
  <div class="mt-2 max-w-full max-h-full overflow-auto">
    <DataTable
      v-model:filters="filters"
      v-model:selection="selection"
      v-model:sortField="sortField"
      :loading="!data.length || props.loading"
      :meta-key-selection="false"
      :selection-mode="props.multiSelect ? 'multiple' : 'single'"
      :value="filteredRows"
      data-key="index"
      filter-display="menu"
      reorderable-columns
      removable-sort
      @column-reorder="onColumnReorder"
      @row-reorder="onRowReorder"
    >
      <template #header>
        <div v-if="activeFilters.length" class="flex">
          <Button
            icon="pi pi-filter-slash"
            label="Clear Filters"
            class="p-button-outlined p-button-sm"
            @click="resetFilters()"
          />
          <Chip
            v-for="filter in activeFilters"
            :key="filter.text"
            :label="filter.text"
            class="ml-2"
            removable
            @remove="resetFilter(filter)"
          />
        </div>
      </template>
      <template #empty><div class="p-4">No entries</div></template>
      <template #loading><div class="p-4">Loading...</div></template>

      <Column
        v-if="props.reorderable"
        :reorderable-column="false"
        header-class="w-3rem"
        row-reorder
      />
      <Column
        v-if="props.multiSelect"
        selection-mode="multiple"
        header-class="w-3rem"
      >
        <template #body="{ data }">
          <Checkbox v-model="selection" :value="data" @click.stop />
        </template>
      </Column>
      <Column field="index" header="Frame" body-class="text-center" sortable>
        <template #body="{ data }">{{ data.index + 1 }}</template>
      </Column>
      <Column field="image" header="Image" body-class="py-0 text-center">
        <template #body="{ data }">
          <Image
            v-if="getImage(data.index)"
            :src="`data:image/png;base64,${getImage(data.index)}`"
            image-class="h-4rem"
            preview
          />
          <div v-else class="h-4rem inline-flex align-items-center">
            <i class="pi pi-exclamation-triangle text-500 text-5xl" />
          </div>
        </template>
      </Column>
      <Column
        v-for="col in filteredColumns"
        :key="col"
        :field="col"
        :header="col"
        :data-type="props.columnTypes[col]"
        body-class="text-center"
        sortable
      >
        <template #filter="{ filterCallback, filterModel }">
          <InputText
            v-if="columnTypes[col] === 'text'"
            v-model="filterModel.value"
            :placeholder="`Search ${col}`"
            class="p-column-filter"
            @keydown.enter="filterCallback"
          />
          <InputNumber
            v-if="columnTypes[col] === 'numeric'"
            v-model="filterModel.value"
            class="p-column-filter"
            @keydown.enter="filterCallback"
          />
        </template>
      </Column>
    </DataTable>
  </div>
</template>

<style scoped>
:deep(.p-datatable-thead) {
  position: sticky;
  top: 0;
  z-index: 1;
}
</style>
