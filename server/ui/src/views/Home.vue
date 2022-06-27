<script setup>
import { onMounted, ref } from 'vue'
import { useSessionStore } from '../store/session'

import IntroPanel from '../components/IntroPanel.vue'

const id = ref('')
const session = useSessionStore()

onMounted(() => {
  session.disconnect()
})
</script>

<template>
  <IntroPanel>
    <form class="flex flex-column" @submit.prevent="$router.push(`/${id}`)">
      <div class="mb-3 text-xl">Enter an ID to connect:</div>
      <InputText v-model="id" class="mb-2" type="text" />
      <Button class="text-center" type="submit">Connect</Button>
    </form>
  </IntroPanel>
</template>
