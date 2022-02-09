import { createApp } from 'vue'
import App from './App.vue'
import PrimeVue from 'primevue/config'
import ToastService from 'primevue/toastservice'
import router from './router'

import Button from 'primevue/button'
import Checkbox from 'primevue/checkbox'
import Dropdown from 'primevue/dropdown'
import InputText from 'primevue/inputtext'
import Menu from 'primevue/menu'
import MultiSelect from 'primevue/multiselect'
import Toast from 'primevue/toast'

import 'primevue/resources/themes/arya-blue/theme.css'
import 'primevue/resources/primevue.min.css'
import 'primeflex/primeflex.css'
import 'primeicons/primeicons.css'

const app = createApp(App)

app.component('Button', Button)
app.component('Checkbox', Checkbox)
app.component('Dropdown', Dropdown)
app.component('InputText', InputText)
app.component('Menu', Menu)
app.component('MultiSelect', MultiSelect)
app.component('Toast', Toast)

app.use(PrimeVue, {
  ripple: true,
  inputStyle: 'filled'
})
app.use(ToastService)
app.use(router)

app.mount('#app')
