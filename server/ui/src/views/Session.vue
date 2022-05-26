<script setup>
import { computed, ref, watch } from 'vue'
import { EVENT, STATUS, useWS } from '../ws'

import ComplexGraph from '../components/ComplexGraph.vue'
import ComplexTable from '../components/ComplexTable.vue'
import ExampleImage from '../assets/image.png'

const props = defineProps({
  id: String
})

const ws = useWS(props.id)
const { complexes, data, images, status } = ws
const loading = ref(true)
const columns = ref([])
const columnTypes = ref({})

const selectedColumns = ref([...columns.value])
const selectedComplex = ref(null)
const selectedFrame = ref(null)
const selectedRows = ref([])
const hiddenFrames = ref([])

const selectionMode = ref(false)
const reorderMode = ref(false)
const nameColumn = ref(null)
const oldData = ref([])

const displayColumns = computed(() => {
  return columns.value.filter(c => selectedColumns.value.includes(c))
})

const selectedIndices = computed(() => {
  return selectedRows.value.map(r => r.index)
})

watch(data, () => {
  // trim string columns
  data.value.forEach(o => {
    Object.entries(o).forEach(([key, value]) => {
      if (typeof value !== 'string') return
      o[key] = value.trim()
    })
  })

  const columnSet = new Set([].concat(...data.value.map(o => Object.keys(o))))
  columnSet.delete('index')
  columns.value = Array.from(columnSet)

  const name = columns.value.find(c => c.toLowerCase() === 'name')
  if (name) nameColumn.value = name
})

watch(
  columns,
  () => {
    selectedColumns.value = [...columns.value]

    // hide columns that have > 30 char data
    data.value.forEach(o => {
      Object.entries(o).forEach(([key, value]) => {
        if (value.length < 30) return
        const index = selectedColumns.value.indexOf(key)
        if (index !== -1) selectedColumns.value.splice(index, 1)
      })
    })

    // guess types of columns based on data
    columnTypes.value = {}
    data.value.forEach(o => {
      Object.entries(o).forEach(([key, value]) => {
        if (columnTypes.value[key] === 'text') return
        const type = isNaN(value) ? 'text' : 'numeric'
        columnTypes.value[key] = type
      })
    })
  },
  { immediate: true }
)

ws.on(EVENT.DATA, () => {
  loading.value = false
})

ws.on(EVENT.SELECT_COMPLEX, value => {
  selectedComplex.value = value
})

ws.on(EVENT.SELECT_FRAME, index => {
  selectedFrame.value = data.value.find(c => c.index === index)
})

const deleteSelection = () => {
  loading.value = true
  ws.send(EVENT.DELETE_FRAMES, selectedIndices.value)
  selectedRows.value = []
}

const hideSelection = () => {
  hiddenFrames.value.push(...selectedIndices.value)
  selectedRows.value = []
}

const saveReorder = () => {
  loading.value = true
  oldData.value = [...data.value]
  const ids = data.value.map(f => f.index)
  ws.send(EVENT.REORDER_FRAMES, ids)
  reorderMode.value = false
}

const selectComplex = event => {
  data.value = []
  selectedFrame.value = null
  if (!event) return
  ws.send(EVENT.SELECT_COMPLEX, event.value)
}

const selectFrame = event => {
  if (!event) return
  ws.send(EVENT.SELECT_FRAME, event.index)
}

const splitSelection = (single, remove) => {
  loading.value = remove
  ws.send(EVENT.SPLIT_FRAMES, {
    indices: selectedIndices.value,
    name_column: nameColumn.value,
    single,
    remove
  })
  selectedRows.value = []
}

const toggleReorderMode = () => {
  reorderMode.value = !reorderMode.value
  if (reorderMode.value) {
    oldData.value = [...data.value]
  } else {
    data.value = [...oldData.value]
  }
}

const toggleSelectionMode = () => {
  selectionMode.value = !selectionMode.value
  selectedRows.value = []
}

const unhideAll = () => {
  hiddenFrames.value = []
}

ws.connect()
</script>

<template>
  <div v-if="status === STATUS.ONLINE && !selectedComplex" class="h-full grid">
    <div class="col-5 flex align-items-center justify-content-center">
      <div class="inline-block">
        <h1 class="mt-0 mb-2 text-6xl">Data Table Plugin</h1>
        <div class="mb-8 text-500 text-xl">
          View your multi-frame molecule metadata<br />
          in an interactive table.
        </div>

        <div
          class="surface-card py-5 border-2 border-100 text-center"
          style="border-radius: 8px"
        >
          <div class="mb-3 text-xl">Select an entry to begin</div>
          <Dropdown
            v-model="selectedComplex"
            :options="complexes"
            class="w-15rem"
            option-label="name"
            option-value="index"
            placeholder="click here"
            @change="selectComplex"
          />
        </div>
      </div>
    </div>

    <div
      class="col-7 px-6 flex align-items-center justify-content-center surface-50 text-center"
    >
      <div class="inline-block">
        <h2 class="text-4xl">Example Output</h2>
        <img :src="ExampleImage" alt="example" class="w-full" />
        <div class="mt-4 text-500 text-xl">
          Best used for analyzing metadata for a multi-frame SDF small molecule
          entry
        </div>
      </div>
    </div>
  </div>

  <div v-else class="h-full flex align-items-center justify-content-center">
    <div
      class="flex flex-column surface-card max-w-full max-h-full min-w-2 p-4 shadow-2 border-round text-center"
    >
      <div class="absolute bottom-0 right-0 p-2 text-300">
        {{ props.id }}
      </div>

      <div v-if="status === STATUS.OFFLINE">
        <div>Disconnected</div>
        <div class="py-4">
          <i class="pi pi-times-circle text-6xl" />
        </div>
        <Button @click="ws.connect">Reconnect</Button>
      </div>

      <div v-else-if="status === STATUS.CONNECTING">
        <div>Connecting...</div>
        <div class="py-4">
          <i class="pi pi-spin pi-spinner text-6xl" />
        </div>
      </div>

      <template v-else>
        <div>
          <div class="mx-2 inline-block">
            <div class="mb-2 text-sm text-left">Entry</div>
            <Dropdown
              v-model="selectedComplex"
              :options="complexes"
              class="w-15rem"
              option-label="name"
              option-value="index"
              placeholder="select an entry"
              @change="selectComplex"
            />
          </div>

          <div class="mx-2 inline-block">
            <div class="mb-2 text-sm text-left">Show Columns</div>
            <MultiSelect
              v-model="selectedColumns"
              :options="columns"
              :max-selected-labels="0.1"
              class="w-15rem"
              placeholder="toggle columns"
              selected-items-label="toggle columns"
            />
          </div>

          <div class="mx-2 inline-block">
            <div class="mb-2 text-sm text-left">Name Column</div>
            <Dropdown
              v-model="nameColumn"
              :options="columns"
              class="w-15rem"
              placeholder="no name column"
            />
          </div>
        </div>

        <div class="flex flex-grow-1 min-h-0">
          <ComplexTable
            class="flex-grow-1"
            :columns="displayColumns"
            :column-types="columnTypes"
            :complex-index="selectedComplex"
            :hidden-frames="hiddenFrames"
            :images="images"
            :loading="loading"
            :multi-select="selectionMode"
            :name-column="nameColumn"
            :reorderable="reorderMode"
            v-model="data"
            v-model:selectedFrame="selectedFrame"
            v-model:selectedRows="selectedRows"
            @update:selectedFrame="selectFrame"
          />

          <div class="mt-4 p-2 w-12">
            <ComplexGraph
              :columns="columns"
              :column-types="columnTypes"
              :complex-index="selectedComplex"
              :data="data"
              :images="images"
              :name-column="nameColumn"
              v-model:selectedFrame="selectedFrame"
              @update:selectedFrame="selectFrame"
            />
          </div>
        </div>

        <div class="pt-4">
          <template v-if="selectionMode">
            <Button
              :disabled="!selectedRows.length"
              class="mx-2 p-button-danger"
              @click="deleteSelection"
            >
              Delete
            </Button>
            <Menu
              ref="duplicate"
              :model="[
                {
                  label: 'single entry',
                  icon: 'pi pi-file',
                  command: () => splitSelection(true, false)
                },
                {
                  label: 'multiple entries',
                  icon: 'pi pi-copy',
                  command: () => splitSelection(false, false)
                }
              ]"
              popup
            />
            <Button
              :disabled="!selectedRows.length"
              class="mx-2"
              @click="e => $refs.duplicate.toggle(e)"
            >
              Duplicate
            </Button>
            <Menu
              ref="split"
              :model="[
                {
                  label: 'single entry',
                  icon: 'pi pi-file',
                  command: () => splitSelection(true, true)
                },
                {
                  label: 'multiple entries',
                  icon: 'pi pi-copy',
                  command: () => splitSelection(false, true)
                }
              ]"
              popup
            />
            <Button
              :disabled="!selectedRows.length"
              class="mx-2"
              @click="e => $refs.split.toggle(e)"
            >
              Split
            </Button>
            <Button
              :disabled="!selectedRows.length"
              class="mx-2"
              @click="hideSelection"
            >
              Hide
            </Button>
          </template>

          <template v-else-if="reorderMode">
            <Button class="mx-2" @click="saveReorder"> Save </Button>
          </template>

          <Button
            v-if="hiddenFrames.length"
            class="mx-2 p-button-outlined"
            @click="unhideAll"
          >
            Unhide All
          </Button>

          <Button
            v-if="!reorderMode"
            class="mx-2 p-button-outlined"
            @click="toggleSelectionMode"
          >
            {{ selectionMode ? 'Cancel' : 'Selection Mode' }}
          </Button>

          <!-- <Button
            v-if="!selectionMode"
            class="mx-2 p-button-outlined"
            @click="toggleReorderMode"
          >
            {{ reorderMode ? 'Cancel' : 'Reorder Mode' }}
          </Button> -->
        </div>
      </template>
    </div>
  </div>
</template>
