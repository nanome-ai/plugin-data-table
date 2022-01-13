<script setup>
import { ref, watch } from 'vue'
import { useVModel } from '../composables.js'

import DataTable from 'primevue/datatable'
import Column from 'primevue/column'
import Image from 'primevue/image'

const props = defineProps({
  columns: Array,
  complexIndex: Number,
  images: Object,
  modelValue: Array,
  selectedFrame: Object
})

const emit = defineEmits(['update:modelValue', 'update:selectedFrame'])

const data = useVModel(props, emit)
const selectedFrame = useVModel(props, emit, 'selectedFrame')

const selectedColumns = ref([...props.columns])
const sortField = ref(null)

watch(
  () => props.columns,
  () => {
    selectedColumns.value = [...props.columns]

    // hide columns that have > 30 char data
    data.value.forEach(o => {
      Object.entries(o).forEach(([key, value]) => {
        if (value.length > 30) {
          const index = selectedColumns.value.indexOf(key)
          if (index !== -1) selectedColumns.value.splice(index, 1)
        }
      })
    })
  },
  { immediate: true }
)

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
</script>

<template>
  <MultiSelect
    v-model="selectedColumns"
    :options="props.columns"
    :max-selected-labels="0.1"
    class="mx-2"
    placeholder="toggle columns"
    selected-items-label="toggle columns"
  />

  <div class="mt-2 max-w-full max-h-full overflow-x-scroll">
    <DataTable
      v-model:selection="selectedFrame"
      v-model:sortField="sortField"
      :loading="!data.length"
      :value="data"
      selection-mode="single"
      reorderable-columns
      removable-sort
      @column-reorder="onColumnReorder"
      @row-reorder="onRowReorder"
    >
      <Column :reorderable-column="false" header-class="w-3rem" row-reorder />
      <Column field="index" header="Index" body-class="text-center" sortable>
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
        v-for="col in props.columns.filter(c => selectedColumns.includes(c))"
        :key="col"
        :field="col"
        :header="col"
        body-class="text-center"
        sortable
      />
    </DataTable>
  </div>
</template>
