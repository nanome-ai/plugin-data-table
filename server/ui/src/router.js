import { createRouter, createWebHistory } from 'vue-router'
import Home from './views/Home.vue'
import Session from './views/Session.vue'

const routes = [
  {
    path: '/',
    name: 'Home',
    component: Home
  },
  {
    path: '/:id',
    name: 'Session',
    component: Session,
    props: true
  }
]

const router = createRouter({
  history: createWebHistory(),
  routes
})

export default router
