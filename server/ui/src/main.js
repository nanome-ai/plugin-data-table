import { createApp, markRaw } from 'vue'
import { createPinia } from 'pinia'
import App from './App.vue'
import PrimeVue from 'primevue/config'
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
import Dropdown from 'primevue/dropdown'
import InputNumber from 'primevue/inputnumber'
import InputText from 'primevue/inputtext'
import Menu from 'primevue/menu'
import MultiSelect from 'primevue/multiselect'
import Skeleton from 'primevue/skeleton'
import Toast from 'primevue/toast'
import ToggleButton from 'primevue/togglebutton'

const app = createApp(App)
const pinia = createPinia()
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
app.component('Dropdown', Dropdown)
app.component('InputNumber', InputNumber)
app.component('InputText', InputText)
app.component('Menu', Menu)
app.component('MultiSelect', MultiSelect)
app.component('Skeleton', Skeleton)
app.component('Toast', Toast)
app.component('ToggleButton', ToggleButton)

app.use(PrimeVue, {
  ripple: true,
  inputStyle: 'filled'
})
app.use(ToastService)
app.use(router)

app.mount('#app')
