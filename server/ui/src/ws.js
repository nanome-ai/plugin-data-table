export const STATUS = {
  CONNECTING: 'connecting',
  OFFLINE: 'offline',
  ONLINE: 'online'
}

export const EVENT = {
  CLOSE: 'close',
  COMPLEXES: 'complexes',
  DELETE_FRAMES: 'delete-frames',
  ERROR: 'error',
  FRAMES: 'frames',
  IMAGE: 'image',
  JOIN: 'join',
  REORDER_FRAMES: 'reorder-frames',
  SELECT_COMPLEX: 'select-complex',
  SELECT_FRAME: 'select-frame',
  SPLIT_FRAMES: 'split-frames',
  STATUS: 'status'
}

export function useWS(id) {
  let ws = null
  const listeners = {}

  const send = (type, data) => {
    console.log('send', type, data)
    ws.send(JSON.stringify({ type, data }))
  }

  const emit = (type, data) => {
    if (!listeners[type]) return
    listeners[type].forEach(fn => fn(data))
  }

  const onOpen = () => {
    send(EVENT.JOIN, id)
    emit(EVENT.STATUS, STATUS.ONLINE)
  }

  const onMessage = e => {
    const msg = JSON.parse(e.data)
    console.log('recv', msg.type)

    if (msg.type === EVENT.CLOSE) {
      ws.close()
    }

    emit(msg.type, msg.data)
  }

  const onClose = () => {
    ws = null
    emit(EVENT.STATUS, STATUS.OFFLINE)
  }

  const connect = () => {
    if (ws) return
    try {
      const protocol = location.protocol === 'https:' ? 'wss:' : 'ws:'
      ws = new WebSocket(`${protocol}//${location.host}/ws`)
      emit(EVENT.STATUS, STATUS.CONNECTING)

      ws.onopen = onOpen
      ws.onmessage = onMessage
      ws.onclose = onClose
    } catch (e) {
      console.error(e)
    }
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
    connect,
    send,
    on
  }
}
