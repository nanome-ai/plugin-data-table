const WebSocket = require('ws')
let wss = null
let sessions = {}

const send = (ws, type, data) => {
  ws.send(JSON.stringify({ type, data }))
}

const broadcast = (ws, type, data, toPlugin = null) => {
  const logType = toPlugin
    ? 'send_plugin'
    : toPlugin !== null
    ? 'send_clients'
    : 'broadcast'

  if (type !== 'image') {
    const debug = type === 'data' ? '---' : data
    console.log(ws.session, logType, type, debug)
  }

  wss.clients.forEach(client => {
    if (client === ws) return
    if (client.session !== ws.session) return
    if (toPlugin !== null && client.plugin !== toPlugin) return
    send(client, type, data)
  })
}

const send_clients = (ws, type, data) => {
  broadcast(ws, type, data, false)
}

const send_plugin = (ws, type, data) => {
  broadcast(ws, type, data, true)
}

const onMessage = ws => e => {
  let msg
  try {
    msg = JSON.parse(e)
  } catch (e) {
    ws.terminate()
    return
  }

  const { type, data } = msg

  if (!ws.session && !['host', 'join'].includes(type)) {
    send(ws, 'error', 'Not connected')
    ws.close()
    return
  }

  switch (type) {
    case 'host':
      ws.plugin = true
      ws.session = data
      sessions[data] = true
      break
    case 'join':
      if (sessions[data]) {
        ws.session = data
        send_plugin(ws, type)
      } else {
        send(ws, 'error', `Invalid session "${data}"`)
        ws.close()
      }
      break
    case 'complexes':
    case 'data':
    case 'image':
      send_clients(ws, type, data)
      break
    case 'delete-frames':
    case 'reorder-frames':
    case 'split-frames':
      send_plugin(ws, type, data)
      break
    case 'select-complex':
    case 'select-frame':
      broadcast(ws, type, data)
      break
    default:
      send(ws, 'error', `Unknown command "${type}"`)
      ws.close()
  }
}

const onClose = ws => () => {
  if (!ws.plugin) return
  delete sessions[ws.session]
  broadcast(ws, 'close')
}

const onConnection = ws => {
  ws.alive = true
  ws.plugin = false
  ws.session = null

  ws.on('message', onMessage(ws))
  ws.on('close', onClose(ws))
  ws.on('pong', () => (ws.alive = true))
}

const killInactive = () => {
  wss.clients.forEach(ws => {
    if (!ws.alive) {
      ws.terminate()
      return
    }
    ws.alive = false
    ws.ping()
  })
}

exports.init = server => {
  wss = new WebSocket.WebSocketServer({ server, path: '/ws' })
  wss.on('connection', onConnection)

  const id = setInterval(killInactive, 30000)
  wss.on('close', () => clearInterval(id))
}
