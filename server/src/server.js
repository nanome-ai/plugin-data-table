const http = require('http')
const app = require('./app')
const ws = require('./ws')

const server = http.createServer(app)
ws.init(server)
server.listen(80)
