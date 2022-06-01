<script setup>
import { ref } from 'vue'
import { STATUS } from '../ws'
import { useSessionStore } from '../store/session'

import ComplexGraph from '../components/ComplexGraph.vue'
import ComplexTable from '../components/ComplexTable.vue'
import IntroPanel from '../components/IntroPanel.vue'

const props = defineProps({
  id: String
})

const session = useSessionStore()

const selectionMode = ref(false)
const reorderMode = ref(false)
const showGraph = ref(false)
// const oldData = ref([])

// const saveReorder = () => {
//   loading.value = true
//   oldData.value = [...data.value]
//   const ids = data.value.map(f => f.index)
//   ws.send(EVENT.REORDER_FRAMES, ids)
//   reorderMode.value = false
// }

// const toggleReorderMode = () => {
//   reorderMode.value = !reorderMode.value
//   if (reorderMode.value) {
//     oldData.value = [...data.value]
//   } else {
//     data.value = [...oldData.value]
//   }
// }

const toggleSelectionMode = () => {
  selectionMode.value = !selectionMode.value
  session.selectedFrames = []
}

session.connect(props.id)
</script>

<template>
  <IntroPanel
    v-if="session.status === STATUS.ONLINE && !session.selectedComplex"
  >
    <div class="mb-3 text-xl">Select an entry to begin</div>
    <Dropdown
      v-model="session.selectedComplex"
      :options="session.complexes"
      class="w-full"
      option-label="name"
      option-value="index"
      placeholder="click here"
      @change="({ value }) => session.selectComplex(value)"
    />
  </IntroPanel>

  <div v-else class="h-full flex align-items-center justify-content-center">
    <div
      class="flex flex-column surface-card max-w-full max-h-full min-w-2 p-4 shadow-2 border-round text-center"
    >
      <div class="absolute bottom-0 right-0 p-2 text-300">
        {{ props.id }}
      </div>

      <div v-if="session.status === STATUS.OFFLINE">
        <div>Disconnected</div>
        <div class="py-4">
          <i class="pi pi-times-circle text-6xl" />
        </div>
        <Button @click="session.connect(props.id)">Reconnect</Button>
      </div>

      <div v-else-if="session.status === STATUS.CONNECTING">
        <div>Connecting...</div>
        <div class="py-4">
          <i class="pi pi-spin pi-spinner text-6xl" />
        </div>
      </div>

      <template v-else>
        <div class="flex flex-grow-1 min-h-0">
          <div class="flex flex-column flex-grow-1 min-w-0">
            <div>
              <div class="mx-2 inline-block">
                <div class="mb-2 text-sm text-left">Entry</div>
                <Dropdown
                  v-model="session.selectedComplex"
                  :options="session.complexes"
                  class="w-15rem"
                  option-label="name"
                  option-value="index"
                  placeholder="select an entry"
                  @change="({ value }) => session.selectComplex(value)"
                />
              </div>

              <div class="mx-2 inline-block">
                <div class="mb-2 text-sm text-left">Show Columns</div>
                <MultiSelect
                  v-model="session.selectedColumns"
                  :options="session.columns"
                  :max-selected-labels="0.1"
                  class="w-15rem"
                  placeholder="toggle columns"
                  selected-items-label="toggle columns"
                />
              </div>

              <div class="mx-2 inline-block">
                <div class="mb-2 text-sm text-left">Name Column</div>
                <Dropdown
                  v-model="session.nameColumn"
                  :options="session.columns"
                  class="w-15rem"
                  placeholder="no name column"
                />
              </div>

              <ToggleButton
                v-model="showGraph"
                class="mx-2"
                on-icon="pi pi-chart-bar"
                on-label="graphs"
                off-icon="pi pi-chart-bar"
                off-label="graphs"
              />
            </div>

            <ComplexTable
              class="flex-grow-1"
              :multi-select="selectionMode"
              :reorderable="reorderMode"
            />
          </div>

          <div v-if="showGraph" class="pl-4 w-12">
            <div class="mb-6 text-xl">Quick Graphs</div>
            <ComplexGraph />
          </div>
        </div>

        <div class="pt-4">
          <template v-if="selectionMode">
            <Button
              :disabled="!session.selectedFrames.length"
              class="mx-2 p-button-danger"
              @click="session.deleteSelection"
            >
              Delete
            </Button>
            <Menu
              ref="duplicate"
              :model="[
                {
                  label: 'single entry',
                  icon: 'pi pi-file',
                  command: () => session.splitSelection(true, false)
                },
                {
                  label: 'multiple entries',
                  icon: 'pi pi-copy',
                  command: () => session.splitSelection(false, false)
                }
              ]"
              popup
            />
            <Button
              :disabled="!session.selectedFrames.length"
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
                  command: () => session.splitSelection(true, true)
                },
                {
                  label: 'multiple entries',
                  icon: 'pi pi-copy',
                  command: () => session.splitSelection(false, true)
                }
              ]"
              popup
            />
            <Button
              :disabled="!session.selectedFrames.length"
              class="mx-2"
              @click="e => $refs.split.toggle(e)"
            >
              Split
            </Button>
            <Button
              :disabled="!session.selectedFrames.length"
              class="mx-2"
              @click="session.hideSelection"
            >
              Hide
            </Button>
          </template>

          <template v-else-if="reorderMode">
            <Button class="mx-2" @click="saveReorder"> Save </Button>
          </template>

          <Button
            v-if="session.hiddenFrames.length"
            class="mx-2 p-button-outlined"
            @click="session.unhideAll"
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
