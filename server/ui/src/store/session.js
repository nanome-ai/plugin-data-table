import { defineStore } from 'pinia'
import { STATUS, EVENT, useWS } from '../ws'

const setupEventListeners = (store, ws) => {
  ws.on(EVENT.STATUS, status => {
    store.status = status
  })

  ws.on(EVENT.COMPLEXES, complexes => {
    store.complexes = complexes
  })

  ws.on(EVENT.FRAMES, frames => {
    // trim string columns
    frames.forEach(o => {
      Object.entries(o).forEach(([key, value]) => {
        if (typeof value !== 'string') return
        o[key] = value.trim()
      })
    })

    // get unique column names
    const columnSet = new Set([].concat(...frames.map(o => Object.keys(o))))
    columnSet.delete('index')
    const columns = Array.from(columnSet)

    // guess types of columns based on data
    const columnTypes = {}
    frames.forEach(o => {
      Object.entries(o).forEach(([key, value]) => {
        if (columnTypes[key] === 'text') return
        const type = isNaN(value) ? 'text' : 'numeric'
        columnTypes[key] = type
      })
    })

    // hide columns that have > 30 char data
    const selectedColumns = [...columns]
    frames.forEach(o => {
      Object.entries(o).forEach(([key, value]) => {
        if (value.length < 30) return
        const index = selectedColumns.indexOf(key)
        if (index !== -1) selectedColumns.splice(index, 1)
      })
    })

    const name = columns.find(c => c.toLowerCase() === 'name')
    if (name) store.nameColumn = name

    store.frames = frames
    store.columns = columns
    store.columnTypes = columnTypes
    store.selectedColumns = selectedColumns
  })

  ws.on(EVENT.IMAGE, ({ id, data }) => {
    store.images[id] = data
  })

  ws.on(EVENT.SELECT_COMPLEX, index => {
    store.selectedComplex = index
  })

  ws.on(EVENT.SELECT_FRAME, index => {
    store.selectedFrame = store.frames.find(f => f.index === index)
  })
}

export const useSessionStore = defineStore('session', {
  state: () => ({
    ws: null,
    status: STATUS.OFFLINE,
    loading: false,
    complexes: [],
    images: {},
    frames: [],
    columns: [],
    columnTypes: {},
    nameColumn: null,
    hiddenFrames: [],
    selectedComplex: null,
    selectedColumns: [],
    selectedFrame: null,
    selectedFrames: []
  }),

  getters: {
    displayColumns() {
      return this.columns.filter(c => this.selectedColumns.includes(c))
    },

    getImage() {
      const complexId = this.selectedComplex
      const prefix = 'data:image/png;base64,'
      return id => {
        const image = this.images[`${complexId}-${id}`]
        if (!image) return null
        return prefix + image
      }
    },

    numericColumns() {
      return this.columns.filter(c => this.columnTypes[c] === 'numeric')
    },

    selectedFrameIndices() {
      return this.selectedFrames.map(f => f.index)
    }
  },

  actions: {
    connect(id) {
      const ws = useWS(id)
      this.ws = ws
      setupEventListeners(this, ws)
      ws.connect()
    },

    disconnect() {
      if (!this.ws) return
      this.ws.disconnect()
    },

    deleteSelection() {
      this.loading = true
      this.ws.send(EVENT.DELETE_FRAMES, this.selectedFrameIndices)
      this.selectedFrames = []
    },

    hideSelection() {
      this.hiddenFrames.push(...this.selectedFrameIndices)
      this.selectedFrames = []
    },

    selectComplex(index) {
      this.frames = []
      this.selectedFrame = null
      this.ws.send(EVENT.SELECT_COMPLEX, index)
    },

    selectFrame(index) {
      this.selectedFrame = this.frames.find(f => f.index === index)
      this.ws.send(EVENT.SELECT_FRAME, index)
    },

    splitSelection(single, remove) {
      this.loading = remove
      this.ws.send(EVENT.SPLIT_FRAMES, {
        indices: this.selectedFrameIndices,
        name_column: this.nameColumn,
        single,
        remove
      })
      this.selectedFrames = []
    },

    unhideAll() {
      this.hiddenFrames = []
    }
  }
})
