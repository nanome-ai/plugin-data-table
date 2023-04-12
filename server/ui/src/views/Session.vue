<script setup>
import { computed, ref } from 'vue'
import { useConfirm } from 'primevue/useconfirm'
import { STATUS } from '../ws'
import { useSessionStore } from '../store/session'

import ComplexGraph from '../components/ComplexGraph.vue'
import ComplexTable from '../components/ComplexTable.vue'
import EditFrame from '../components/EditFrame.vue'
import Ketcher from '../components/Ketcher.vue'
import NewColumn from '../components/NewColumn.vue'
import IntroPanel from '../components/IntroPanel.vue'

const props = defineProps({
  id: String
})

const session = useSessionStore()
const confirm = useConfirm()

const settings = ref(null)
const showGraphs = ref(false)

const showFullscreenGraph = computed({
  get: () => session.selectedGraph !== null,
  set: value => {
    if (!value) session.selectGraph(null)
  }
})

const toggleSelectionMode = () => {
  session.selectionMode = !session.selectionMode
  session.selectedFrames = []
}

const confirmDelete = e => {
  confirm.require({
    target: e.currentTarget,
    message: `Are you sure you want to delete the selected frame(s)? \nIf there are no frames left, the entry will be removed from Nanome.`,
    icon: 'pi pi-trash',
    acceptClass: 'p-button-danger',
    acceptLabel: 'Delete',
    rejectClass: 'p-button-secondary p-button-text',
    rejectLabel: 'Cancel',
    accept: () => {
      session.deleteSelection()
    }
  })
}

session.connect(props.id)
</script>

<template>
  <IntroPanel
    v-if="session.status === STATUS.ONLINE && !session.selectedComplexes.length"
  >
    <div class="mb-3 text-xl">Select an entry to begin</div>
    <MultiSelect
      v-model="session.selectedComplexes"
      :options="session.complexes"
      class="w-full"
      option-label="name"
      option-value="index"
      placeholder="click here"
      @change="e => session.selectComplexes(e.value)"
    />
  </IntroPanel>

  <div v-else class="h-full flex-center">
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
        <div class="flex min-h-0">
          <div class="flex flex-column min-w-0">
            <div class="mt-2 flex justify-content-center gap-2">
              <span class="p-float-label">
                <MultiSelect
                  v-model="session.selectedComplexes"
                  :options="session.complexes"
                  class="w-15rem p-inputwrapper-filled"
                  option-label="name"
                  option-value="index"
                  placeholder="select an entry"
                  @change="e => session.selectComplexes(e.value)"
                />
                <label>Entries</label>
              </span>

              <NewColumn />

              <Button
                class="p-button-outlined"
                icon="pi pi-cog"
                label="settings"
                @click="e => settings.toggle(e)"
              />

              <OverlayPanel ref="settings">
                <div class="mt-1 flex flex-column gap-3">
                  <span class="mt-3 p-float-label">
                    <MultiSelect
                      :model-value="session.selectedColumns"
                      :options="session.columns"
                      :max-selected-labels="0.1"
                      class="w-full p-inputwrapper-filled"
                      placeholder="toggle columns"
                      selected-items-label="toggle columns"
                      @change="e => session.selectColumns(e.value)"
                    />
                    <label>Show Columns</label>
                  </span>

                  <span class="mt-3 p-float-label">
                    <Dropdown
                      v-model="session.nameColumn"
                      :options="session.columns"
                      class="w-full p-inputwrapper-filled"
                      placeholder="no name column"
                    />
                    <label>Name Column</label>
                  </span>

                  <label class="flex align-items-center">
                    <i class="mr-2 pi pi-server" />
                    <span class="mr-auto">Auto Calc Properties</span>
                    <InputSwitch v-model="session.autoCalcProperties" />
                  </label>

                  <label class="flex align-items-center">
                    <i class="mr-2 pi pi-eye" />
                    <span class="mr-auto">Large Thumbnails</span>
                    <InputSwitch v-model="session.largeThumbnails" />
                  </label>

                  <label class="flex align-items-center">
                    <i class="mr-2 pi pi-search-plus" />
                    <span class="mr-auto">Font Size</span>
                    <SelectButton
                      class="ml-3"
                      v-model="session.fontSize"
                      :options="[
                        { label: 'sm', value: '12px' },
                        { label: 'md', value: '16px' },
                        { label: 'lg', value: '24px' }
                      ]"
                      option-label="label"
                      option-value="value"
                    />
                  </label>

                  <Button
                    v-if="!session.autoCalcProperties"
                    v-tooltip.bottom="'using RDKit'"
                    :disabled="session.hasRDKitProperties"
                    :loading="session.loading"
                    class="p-button-outlined"
                    label="calculate properties"
                    icon="pi pi-server"
                    @click="session.calculateProperties"
                  />
                </div>
              </OverlayPanel>

              <ToggleButton
                v-model="showGraphs"
                class="ml-auto"
                on-icon="pi pi-chart-bar"
                off-icon="pi pi-chart-bar"
                on-label="hide graphs"
                off-label="show graphs"
              />
            </div>

            <ComplexTable
              class="flex-grow-1"
              :multi-select="session.selectionMode"
            />
          </div>

          <template v-if="showGraphs">
            <Divider layout="vertical" />

            <div
              class="pl-2 flex flex-column overflow-auto"
              style="min-width: 33%"
            >
              <ComplexGraph
                v-for="graph in session.graphs"
                :key="graph.id"
                :graph="graph"
              />

              <div class="flex-center flex-wrap gap-2">
                <Button
                  class="p-button-text"
                  label="new graph"
                  icon="pi pi-plus"
                  @click="session.addGraph(false)"
                />

                <Button
                  v-if="session.selectedFrames.length"
                  class="p-button-text"
                  label="new graph from selection"
                  icon="pi pi-plus"
                  @click="session.addGraph(true)"
                />
              </div>
            </div>
          </template>
        </div>

        <div class="pt-4">
          <template v-if="session.selectionMode">
            <Button
              :disabled="!session.selectedFrames.length"
              class="mx-2 p-button-danger"
              icon="pi pi-trash"
              label="Delete"
              @click="confirmDelete"
            />
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
              icon="pi pi-clone"
              label="Duplicate"
              @click="e => $refs.duplicate.toggle(e)"
            />
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
              icon="pi pi-clone"
              label="Split"
              @click="e => $refs.split.toggle(e)"
            />
            <Button
              :disabled="!session.selectedFrames.length"
              class="mx-2"
              icon="pi pi-eye-slash"
              label="Hide"
              @click="session.hideSelection"
            />
          </template>

          <Button
            v-if="session.hiddenFrames.length"
            class="mx-2 p-button-outlined"
            icon="pi pi-eye"
            label="Unhide All"
            @click="session.unhideAll"
          />

          <Button
            class="mx-2 p-button-outlined"
            :icon="session.selectionMode ? 'pi pi-times' : 'pi pi-check-square'"
            :label="session.selectionMode ? 'Cancel' : 'Selection Mode'"
            @click="toggleSelectionMode"
          />

          <template v-if="!session.selectionMode">
            <EditFrame />
            <Ketcher />
          </template>
        </div>

        <Sidebar v-model:visible="showFullscreenGraph" position="full">
          <ComplexGraph
            v-if="showFullscreenGraph"
            :graph="session.selectedGraph"
            fullscreen
          />
        </Sidebar>
      </template>
    </div>
  </div>
</template>
