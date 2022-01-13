const express = require('express')
const ospath = require('path')
const app = express()

const STATIC_DIR = ospath.resolve('ui/dist')

app.use(require('morgan')('dev'))
app.use(require('helmet')({ contentSecurityPolicy: false }))

app.use(express.static(STATIC_DIR))
app.get('*', (req, res) => {
  res.sendFile(STATIC_DIR + '/index.html')
})

app.use((error, req, res, next) => {
  res.status(error.status || 500)
  res.json({ error: error.message })
})

module.exports = app
