import { markRaw } from 'vue'
import { defineStore } from 'pinia'
import { STATUS, EVENT, useWS } from '../ws'
import { randStr } from '../utils'

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
        if (type === 'numeric') {
          o[key] = +value
        }
      })
    })

    if (!store.selectedColumns.length) {
      const selectedColumns = [...columns]
      // hide columns that have > 30 char data
      frames.forEach(o => {
        Object.entries(o).forEach(([key, value]) => {
          if (columnTypes[key] === 'numeric' || value.length < 30) return
          const index = selectedColumns.indexOf(key)
          if (index !== -1) selectedColumns.splice(index, 1)
        })
      })
      store.selectedColumns = selectedColumns
    }

    const name = columns.find(c => c.toLowerCase() === 'name')
    if (name) store.nameColumn = name

    store.frames = frames
    store.columns = columns
    store.columnTypes = columnTypes
    store.loading = false
  })

  ws.on(EVENT.IMAGE, ({ id, data }) => {
    store.images[id] = data
  })

  ws.on(EVENT.SELECT_COMPLEX, index => {
    store.selectedComplex = index
    store.selectedColumns = []
  })

  ws.on(EVENT.SELECT_FRAME, index => {
    store.selectFrame(index, false)
  })

  ws.on(EVENT.UPDATE_FRAME, data => {
    store.updateFrame(data, false)
  })

  ws.on(EVENT.ERROR, error => {
    store.$toast.add({
      severity: 'error',
      summary: 'Error',
      detail: error,
      life: 5000
    })
    store.$router.push('/')
  })
}

const createGraph = () => ({
  id: randStr(6),
  xColumn: null,
  yColumn: null,
  rColumns: [],
  frames: [],
  type: 'scatter',
  reg: {
    type: 'none',
    order: 2,
    color: '0080ff',
    error: false,
    eq: '',
    r2: 0
  }
})

export const useSessionStore = defineStore('session', {
  state: () => ({
    columns: [],
    columnTypes: {},
    complexes: [],
    frames: [],
    graphs: [createGraph()],
    hiddenFrames: [],
    images: {},
    largeThumbnails: false,
    loading: false,
    nameColumn: null,
    selectedColumns: [],
    selectedComplex: null,
    selectedFrame: null,
    selectedFrames: [],
    selectedGraph: null,
    selectionMode: false,
    status: STATUS.OFFLINE,
    ws: null
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

    hasRDKitProperties() {
      const columns = ['MW', 'logP', 'TPSA', 'HBA', 'HBD', 'RB', 'AR']
      return columns.every(c => this.columns.includes(c))
    },

    numericColumns() {
      return this.columns.filter(c => this.columnTypes[c] === 'numeric')
    },

    selectedFrameIndices() {
      return this.selectedFrames.map(f => f.index)
    }
  },

  actions: {
    // #region websockets
    connect(id) {
      const ws = useWS(id)
      this.ws = markRaw(ws)
      setupEventListeners(this, ws)
      ws.connect()
    },

    disconnect() {
      if (!this.ws) return
      this.ws.disconnect()
    },
    // #endregion

    // #region columns
    addColumn(name, values) {
      this.selectedColumns.push(name)
      this.ws.send(EVENT.ADD_COLUMN, { name, values })
    },

    calculateProperties() {
      this.selectedColumns.push('MW', 'logP', 'TPSA', 'HBA', 'HBD', 'RB', 'AR')
      this.loading = true
      this.ws.send(EVENT.CALCULATE_PROPERTIES)
    },
    // #endregion

    // #region graphs
    addGraph(selectedOnly) {
      const graph = createGraph()
      if (selectedOnly) {
        graph.frames = this.selectedFrameIndices
      }
      this.graphs.push(graph)
    },

    removeGraph(graph) {
      const index = this.graphs.indexOf(graph)
      if (index === -1) return
      this.graphs.splice(index, 1)
      this.selectedGraph = null
    },

    selectGraph(graph) {
      this.selectedGraph = graph
    },
    // #endregion

    // #region selection
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

    selectFrame(index, send = true) {
      this.selectedFrame = this.frames.find(f => f.index === index)
      if (send) this.ws.send(EVENT.SELECT_FRAME, index)
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

    updateFrame(data, send = true) {
      const frame = this.frames.find(f => f.index === data.index)
      if (!frame) return
      Object.assign(frame, data)
      if (send) this.ws.send(EVENT.UPDATE_FRAME, data)
    },

    unhideAll() {
      this.hiddenFrames = []
    }
    // #endregion
  }
})
