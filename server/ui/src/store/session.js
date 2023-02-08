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
    ;['id', 'index', 'frame'].forEach(c => columnSet.delete(c))
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
          const hide =
            key.length > 10 ||
            columnTypes[key] !== 'numeric' ||
            value.length > 30
          if (!hide) return

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

  ws.on(EVENT.SELECT_COMPLEXES, indices => {
    store.selectedComplexes = indices
    store.selectedColumns = []
  })

  ws.on(EVENT.SELECT_FRAME, id => {
    store.selectFrame(id, false)
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
    fontSize: '16px',
    frames: [],
    graphs: [createGraph()],
    hiddenFrames: [],
    images: {},
    largeThumbnails: false,
    loading: false,
    nameColumn: null,
    selectedColumns: [],
    selectedComplexes: [],
    selectedFrame: null,
    selectedFrames: [],
    selectedGraph: null,
    selectionMode: false,
    status: STATUS.OFFLINE,
    ws: null
  }),

  persist: {
    enabled: true,
    strategies: [
      {
        key: 'data-table-settings',
        storage: localStorage,
        paths: ['fontSize', 'largeThumbnails']
      }
    ]
  },

  getters: {
    displayColumns() {
      return this.columns.filter(c => this.selectedColumns.includes(c))
    },

    getImage() {
      const prefix = 'data:image/png;base64,'
      return id => {
        const image = this.images[id]
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

    selectedFrameIds() {
      return this.selectedFrames.map(f => f.id)
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
        graph.frames = this.selectedFrameIds
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
      this.ws.send(EVENT.DELETE_FRAMES, this.selectedFrameIds)
      this.selectedFrames = []
    },

    hideSelection() {
      this.hiddenFrames.push(...this.selectedFrameIds)
      this.selectedFrames = []
    },

    selectComplexes(indices) {
      this.frames = []
      this.selectedFrame = null
      this.ws.send(EVENT.SELECT_COMPLEXES, indices)
    },

    selectFrame(id, send = true) {
      this.selectedFrame = this.frames.find(f => f.id === id)
      if (send) this.ws.send(EVENT.SELECT_FRAME, id)
    },

    splitSelection(single, remove) {
      this.loading = remove
      this.ws.send(EVENT.SPLIT_FRAMES, {
        ids: this.selectedFrameIds,
        name_column: this.nameColumn,
        single,
        remove
      })
      this.selectedFrames = []
    },

    updateFrame(data, send = true) {
      const frame = this.frames.find(f => f.id === data.id)
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
