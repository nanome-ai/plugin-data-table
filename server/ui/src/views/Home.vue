<script setup>
import { ref } from 'vue'
const id = ref('')
</script>

<template>
  <div class="h-full flex align-items-center justify-content-center">
    <div class="surface-card w-20rem p-4 shadow-sm border-round">
      <form class="flex flex-column" @submit.prevent="$router.push(`/${id}`)">
        <div class="mb-2">Enter an ID to connect:</div>
        <InputText v-model="id" class="mb-2" type="text" />
        <Button class="text-center" type="submit">Connect</Button>
      </form>
    </div>
  </div>
</template>
