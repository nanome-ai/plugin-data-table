const fs = require('fs')
const http = require('http')
const https = require('https')

const app = require('./app')
const ws = require('./ws')
ws.init()


const http_port = process.env.HTTP_PORT || 80
const https_port = process.env.HTTPS_PORT || 443

const httpServer = http.createServer(app)
httpServer.on('upgrade', ws.onUpgrade)
httpServer.listen(http_port)

const options = {
  key: fs.readFileSync('./certs/local.key'),
  cert: fs.readFileSync('./certs/local.crt')
}
const httpsServer = https.createServer(options, app)
httpsServer.on('upgrade', ws.onUpgrade)
httpsServer.listen(https_port)
