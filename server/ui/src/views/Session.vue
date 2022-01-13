<script setup>
import { reactive, ref, watch } from 'vue'
import { useRouter } from 'vue-router'
import { useToast } from 'primevue/usetoast'

import ComplexTable from '../components/ComplexTable.vue'

const props = defineProps({
  id: String
})

const toast = useToast()
const router = useRouter()

let ws = null
const state = ref('offline')
const data = ref([])
const images = reactive({})
const columns = ref([])

const complexes = ref([])
const selectedComplex = ref(null)
const selectedFrame = ref(null)

watch(data, () => {
  const columnSet = new Set([].concat(...data.value.map(o => Object.keys(o))))
  columnSet.delete('index')
  columns.value = Array.from(columnSet)
})

const send = (type, data) => {
  ws.send(JSON.stringify({ type, data }))
}

const onWSOpen = () => {
  send('join', props.id)
  state.value = 'online'
}

const onWSMessage = e => {
  const msg = JSON.parse(e.data)
  switch (msg.type) {
    case 'complexes':
      complexes.value = msg.data
      break
    case 'data':
      data.value = msg.data
      break
    case 'image':
      images[msg.data.id] = msg.data.data
      break
    case 'select-complex':
      selectedComplex.value = msg.data
      break
    case 'select-frame':
      selectedFrame.value = data.value.find(c => c.index === msg.data)
      break
    case 'close':
      ws.close()
      break
    case 'error':
      toast.add({
        severity: 'error',
        summary: 'Error',
        detail: msg.data,
        life: 5000
      })
      router.push('/')
      break
    default:
      console.error('Unknown message', msg)
  }
}

const onWSClose = () => {
  state.value = 'offline'
  ws = null
  data.value = []
  selectedFrame.value = null
}

const connect = () => {
  if (ws) return
  ws = new WebSocket(`ws://${location.hostname}/ws`)
  state.value = 'connecting'

  ws.onopen = onWSOpen
  ws.onmessage = onWSMessage
  ws.onclose = onWSClose
}

const onSelectComplex = event => {
  data.value = []
  selectedFrame.value = null
  send('select-complex', event.value)
}

const onSelectFrame = event => {
  send('select-frame', event.index)
}

connect()
</script>

<template>
  <div class="h-full p-4 flex align-items-center justify-content-center">
    <div
      class="surface-card max-w-full max-h-full min-w-2 p-4 shadow-2 border-round text-center overflow-y-auto"
    >
      <template v-if="state === 'offline'">
        <div>Disconnected</div>
        <div class="py-4">
          <i class="pi pi-times-circle text-6xl" />
        </div>
        <Button @click="connect">Reconnect</Button>
      </template>

      <template v-else-if="state === 'connecting'">
        <div>Connecting...</div>
        <div class="py-4">
          <i class="pi pi-spin pi-spinner text-6xl" />
        </div>
      </template>

      <template v-else>
        <Dropdown
          v-model="selectedComplex"
          :options="complexes"
          class="mx-2"
          optionLabel="name"
          optionValue="index"
          placeholder="Select a complex"
          @change="onSelectComplex"
        />
        <ComplexTable
          v-if="selectedComplex"
          :columns="columns"
          :complex-index="selectedComplex"
          :images="images"
          v-model="data"
          v-model:selectedFrame="selectedFrame"
          @update:selectedFrame="onSelectFrame"
        />
      </template>
    </div>
  </div>
</template>
