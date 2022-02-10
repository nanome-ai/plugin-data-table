<script setup>
import { computed, ref } from 'vue'
import { useMux, useVModel } from '../composables.js'

import DataTable from 'primevue/datatable'
import Column from 'primevue/column'
import Image from 'primevue/image'

const props = defineProps({
  columns: Array,
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
  return props.columns.filter(c => c !== props.nameColumn)
})

const filteredRows = computed(() => {
  return data.value.filter(row => !props.hiddenFrames.includes(row.index))
})

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
  <div class="mt-2 max-w-full max-h-full overflow-auto">
    <DataTable
      v-model:selection="selection"
      v-model:sortField="sortField"
      :loading="!data.length || props.loading"
      :value="filteredRows"
      :selection-mode="props.multiSelect ? 'multiple' : 'single'"
      :meta-key-selection="false"
      data-key="index"
      reorderable-columns
      removable-sort
      @column-reorder="onColumnReorder"
      @row-reorder="onRowReorder"
    >
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
        :key="nameColumn"
        :field="nameColumn"
        :header="nameColumn"
        body-class="text-center"
        sortable
      />
      <Column
        v-for="col in filteredColumns"
        :key="col"
        :field="col"
        :header="col"
        body-class="text-center"
        sortable
      />
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
