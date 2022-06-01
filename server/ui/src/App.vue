<script setup>
import { onMounted } from 'vue'

onMounted(() => {
  document.body.addEventListener('click', e => {
    document.body.classList.add('using-mouse')
    if (e.target.tagName === 'BUTTON') e.target.blur()
  })

  document.body.addEventListener('keydown', e => {
    if (e.keyCode === 9) {
      document.body.classList.remove('using-mouse')
    }
  })
})
</script>

<template>
  <router-view />
  <Toast />
</template>

<style>
html {
  background: var(--surface-b);
  color: var(--text-color);
  font-family: var(--font-family);
  font-weight: 400;

  scrollbar-color: var(--surface-b) transparent;
  scrollbar-width: thin;
}

::-webkit-scrollbar {
  width: 10px;
  height: 10px;
  background: transparent;
}

::-webkit-scrollbar-corner {
  background: transparent;
}

::-webkit-scrollbar-thumb {
  background: var(--surface-b);
  border-radius: 5px;
}

body {
  margin: 0;
  overflow-x: hidden;
  overflow-y: auto;
}

body.using-mouse :focus {
  box-shadow: none !important;
  outline: none !important;
}

#app {
  width: 100vw;
  height: 100vh;
  background: var(--surface-b);
  color: var(--text-color);
}
</style>
