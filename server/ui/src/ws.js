import { reactive, ref } from 'vue'
import { useRouter } from 'vue-router'
import { useToast } from 'primevue/usetoast'

export const STATUS = {
  CONNECTING: 'connecting',
  OFFLINE: 'offline',
  ONLINE: 'online'
}

export const EVENT = {
  CLOSE: 'close',
  COMPLEXES: 'complexes',
  DATA: 'data',
  DELETE_FRAMES: 'delete-frames',
  ERROR: 'error',
  IMAGE: 'image',
  JOIN: 'join',
  REORDER_FRAMES: 'reorder-frames',
  SELECT_COMPLEX: 'select-complex',
  SELECT_FRAME: 'select-frame',
  SPLIT_FRAMES: 'split-frames'
}

export function useWS(id) {
  const toast = useToast()
  const router = useRouter()

  let ws = null
  const listeners = {}
  const status = ref(STATUS.OFFLINE)
  const data = ref([])
  const images = reactive({})
  const complexes = ref([])

  const send = (type, data) => {
    console.log('send', type, data)
    ws.send(JSON.stringify({ type, data }))
  }

  const onOpen = () => {
    send(EVENT.JOIN, id)
    status.value = STATUS.ONLINE
  }

  const onMessage = e => {
    const msg = JSON.parse(e.data)
    console.log('recv', msg.type)

    switch (msg.type) {
      case EVENT.CLOSE:
        ws.close()
        break
      case EVENT.ERROR:
        toast.add({
          severity: 'error',
          summary: 'Error',
          detail: msg.data,
          life: 5000
        })
        router.push('/')
        break
      case EVENT.COMPLEXES:
        complexes.value = msg.data
        break
      case EVENT.DATA:
        data.value = msg.data
        break
      case EVENT.IMAGE:
        images[msg.data.id] = msg.data.data
        break
    }

    if (listeners[msg.type]) {
      listeners[msg.type].forEach(fn => fn(msg.data))
    }
  }

  const onClose = () => {
    ws = null
    status.value = STATUS.OFFLINE
  }

  const connect = () => {
    if (ws) return
    try {
      const protocol = location.protocol === 'https:' ? 'wss:' : 'ws:'
      ws = new WebSocket(`${protocol}//${location.host}/ws`)
      status.value = STATUS.CONNECTING

      ws.onopen = onOpen
      ws.onmessage = onMessage
      ws.onclose = onClose
    } catch (e) {}
  }

  const on = (event, handler) => {
    if (!Object.values(EVENT).includes(event)) {
      console.error('Unknown event', event)
      return
    }

    if (!listeners[event]) {
      listeners[event] = []
    }

    listeners[event].push(handler)
  }

  return {
    status,
    data,
    images,
    complexes,
    connect,
    send,
    on
  }
}
