<script setup>
import { computed, ref, watch } from 'vue'
import { useSessionStore } from '../store/session'

import { FilterMatchMode, FilterOperator } from 'primevue/api'
import DataTable from 'primevue/datatable'
import Column from 'primevue/column'

const OPERATOR_MAP = {
  [FilterMatchMode.EQUALS]: '=',
  [FilterMatchMode.NOT_EQUALS]: '≠',
  [FilterMatchMode.LESS_THAN]: '<',
  [FilterMatchMode.LESS_THAN_OR_EQUAL_TO]: '≤',
  [FilterMatchMode.GREATER_THAN]: '>',
  [FilterMatchMode.GREATER_THAN_OR_EQUAL_TO]: '≥'
}

const props = defineProps({
  multiSelect: Boolean
})

const session = useSessionStore()

const sortField = ref(null)
const selection = computed({
  get: () => {
    return props.multiSelect ? session.selectedFrames : session.selectedFrame
  },
  set: value => {
    if (props.multiSelect) {
      session.selectedFrames = value
    } else {
      session.selectedFrame = value
    }
  }
})

const filteredColumns = computed(() => {
  const columns = session.displayColumns.filter(c => c !== session.nameColumn)
  columns.unshift(session.nameColumn)
  return columns
})

const filteredRows = computed(() => {
  return session.frames.filter(f => !session.hiddenFrames.includes(f.id))
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
  for (const column of session.columns) {
    const matchMode = {
      text: FilterMatchMode.CONTAINS,
      numeric: FilterMatchMode.LESS_THAN
    }[session.columnTypes[column]]

    filters.value[column] = {
      operator: FilterOperator.AND,
      constraints: [{ matchMode, value: null }]
    }
  }
}

watch(() => session.columnTypes, resetFilters, { immediate: true })

const onColumnReorder = (...args) => {
  console.log('onColumnReorder', ...args)
}

const onRowSelect = e => {
  if (props.multiSelect) return
  session.selectFrame(e.data.id)
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
  <div class="mt-2 overflow-auto">
    <DataTable
      v-model:filters="filters"
      v-model:selection="selection"
      v-model:sortField="sortField"
      :loading="!session.frames.length || session.loading"
      :meta-key-selection="false"
      :selection-mode="props.multiSelect ? 'multiple' : 'single'"
      :value="filteredRows"
      data-key="id"
      filter-display="menu"
      reorderable-columns
      removable-sort
      @column-reorder="onColumnReorder"
      @row-select="onRowSelect"
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
        v-if="props.multiSelect"
        selection-mode="multiple"
        header-class="w-3rem"
      >
        <template #body="{ data }">
          <Checkbox v-model="selection" :value="data" @click.stop />
        </template>
      </Column>
      <Column
        field="frame"
        header="Frame"
        body-class="text-center white-space-nowrap"
        sortable
      >
        <template #body="{ data }">{{ data.frame }}</template>
      </Column>
      <Column field="image" header="Image" body-class="py-0 text-center">
        <template #body="{ data }">
          <Image
            v-if="session.getImage(data.id)"
            :src="session.getImage(data.id)"
            :image-class="session.largeThumbnails ? 'h-8rem' : 'h-4rem'"
            preview
          />
          <div
            v-else
            :class="session.largeThumbnails ? 'h-8rem' : 'h-4rem'"
            class="inline-flex align-items-center"
          >
            <i class="pi pi-exclamation-triangle text-500 text-5xl" />
          </div>
        </template>
      </Column>
      <Column
        v-for="col in filteredColumns"
        :key="col"
        :field="col"
        :header="col"
        :data-type="session.columnTypes[col]"
        body-class="text-center"
        sortable
      >
        <template #filter="{ filterCallback, filterModel }">
          <InputText
            v-if="session.columnTypes[col] === 'text'"
            v-model="filterModel.value"
            :placeholder="`Search ${col}`"
            class="p-column-filter"
            @keydown.enter="filterCallback"
          />
          <InputNumber
            v-if="session.columnTypes[col] === 'numeric'"
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

:deep(.p-datatable-tbody) {
  user-select: none;
}
</style>
