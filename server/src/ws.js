const WebSocket = require('ws')
let wss = null
let sessions = {}

const send = (ws, type, data) => {
  ws.send(JSON.stringify({ type, data }))
}

const broadcast = (session, type, data, toPlugin = null) => {
  wss.clients.forEach(client => {
    if (client.session !== session) return
    if (toPlugin !== null && client.plugin !== toPlugin) return
    send(client, type, data)
  })
}

const send_clients = (session, type, data) => {
  broadcast(session, type, data, false)
}

const send_plugin = (session, type, data) => {
  broadcast(session, type, data, true)
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
        send_plugin(ws.session, type)
      } else {
        send(ws, 'error', `Invalid session "${data}"`)
        ws.close()
      }
      break
    case 'complexes':
    case 'data':
    case 'image':
      send_clients(ws.session, type, data)
      break
    case 'select-complex':
    case 'select-frame':
      if (ws.plugin) {
        send_clients(ws.session, type, data)
      } else {
        broadcast(ws.session, type, data)
      }
      break
    default:
      send(ws, 'error', `Unknown command "${type}"`)
      ws.close()
  }
}

const onClose = ws => () => {
  if (ws.plugin) {
    delete sessions[ws.session]
    broadcast(ws.session, 'close')
  }
}

const onConnection = ws => {
  ws.alive = true
  ws.plugin = false
  ws.session = null

  ws.on('message', onMessage(ws))
  ws.on('close', onClose(ws))
  ws.on('pong', () => (ws.alive = true))
}

exports.init = server => {
  wss = new WebSocket.WebSocketServer({ server, path: '/ws' })
  wss.on('connection', onConnection)
}
