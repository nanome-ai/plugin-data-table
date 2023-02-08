import { createApp, markRaw } from 'vue'
import { createPinia } from 'pinia'
import piniaPersist from 'pinia-plugin-persist'

import App from './App.vue'
import PrimeVue from 'primevue/config'
import ConfirmationService from 'primevue/confirmationservice'
import ToastEventBus from 'primevue/toasteventbus'
import ToastService from 'primevue/toastservice'
import router from './router'

import 'primevue/resources/themes/arya-blue/theme.css'
import 'primevue/resources/primevue.min.css'
import 'primeflex/primeflex.css'
import 'primeicons/primeicons.css'

import Button from 'primevue/button'
import Checkbox from 'primevue/checkbox'
import Chip from 'primevue/chip'
import ColorPicker from 'primevue/colorpicker'
import Dialog from 'primevue/dialog'
import Divider from 'primevue/divider'
import Dropdown from 'primevue/dropdown'
import Image from 'primevue/image'
import InputNumber from 'primevue/inputnumber'
import InputSwitch from 'primevue/inputswitch'
import InputText from 'primevue/inputtext'
import Menu from 'primevue/menu'
import MultiSelect from 'primevue/multiselect'
import OverlayPanel from 'primevue/overlaypanel'
import SelectButton from 'primevue/selectbutton'
import Sidebar from 'primevue/sidebar'
import Skeleton from 'primevue/skeleton'
import ToggleButton from 'primevue/togglebutton'
import Tooltip from 'primevue/tooltip'

import { Chart } from 'chart.js'
Chart.defaults.color = '#ccc'

const app = createApp(App)
const pinia = createPinia()
pinia.use(piniaPersist)
pinia.use(({ store }) => {
  store.$router = markRaw(router)
  store.$toast = markRaw({
    add: message => ToastEventBus.emit('add', message)
  })
})
app.use(pinia)

app.component('Button', Button)
app.component('Checkbox', Checkbox)
app.component('Chip', Chip)
app.component('ColorPicker', ColorPicker)
app.component('Dialog', Dialog)
app.component('Divider', Divider)
app.component('Dropdown', Dropdown)
app.component('Image', Image)
app.component('InputNumber', InputNumber)
app.component('InputSwitch', InputSwitch)
app.component('InputText', InputText)
app.component('Menu', Menu)
app.component('MultiSelect', MultiSelect)
app.component('OverlayPanel', OverlayPanel)
app.component('SelectButton', SelectButton)
app.component('Sidebar', Sidebar)
app.component('Skeleton', Skeleton)
app.component('ToggleButton', ToggleButton)

app.directive('tooltip', Tooltip)

app.use(PrimeVue, {
  ripple: true
})
app.use(ConfirmationService)
app.use(ToastService)
app.use(router)

app.mount('#app')
