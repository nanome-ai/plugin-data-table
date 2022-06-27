<script setup>
import { computed, ref, watch } from 'vue'
import { useSessionStore } from '../store/session'

const session = useSessionStore()
const visible = ref(false)

const data = ref(null)

const frameJSON = computed(() => {
  return JSON.stringify(session.selectedFrame)
})

const columns = computed(() => {
  const columns = session.columns.filter(c => c !== session.nameColumn)
  columns.unshift(session.nameColumn)
  return columns
})

const hasChanges = computed(() => {
  return frameJSON.value !== JSON.stringify(data.value)
})

const sync = () => {
  data.value = JSON.parse(frameJSON.value)
}

watch(() => session.selectedFrame, sync, { immediate: true })

const selectIndex = dir => {
  const index = session.frames.indexOf(session.selectedFrame)
  let newIndex = index + dir
  if (newIndex < 0) newIndex += session.frames.length
  newIndex %= session.frames.length
  session.selectFrame(newIndex)
}

const save = () => {
  session.updateFrame(data.value)
  visible.value = false
}

const cancel = () => {
  sync()
  visible.value = false
}
</script>

<template>
  <Button class="mx-2 p-button-outlined" @click="visible = true">
    {{ 'Edit Frame' }}
  </Button>

  <Dialog
    v-model:visible="visible"
    style="width: 50vw"
    class="text-center"
    dismissable-mask
    maximizable
    modal
  >
    <template #header>Edit Frame {{ data.index + 1 }}</template>

    <Image
      v-if="session.getImage(data.index)"
      :src="session.getImage(data.index)"
    />
    <div v-else class="h-6rem">
      <i class="pi pi-exclamation-triangle text-500 text-6xl" />
      <div>no image</div>
    </div>

    <div class="mt-5 flex flex-center flex-wrap gap-5">
      <span v-for="column in columns" class="w-15rem p-float-label">
        <InputText
          v-model="data[column]"
          class="w-full p-inputwrapper-filled"
        />
        <label>{{ column }}</label>
      </span>
    </div>

    <template #footer>
      <div class="mt-4 flex">
        <Button
          class="p-button-text"
          icon="pi pi-arrow-left"
          label="prev"
          @click="selectIndex(-1)"
        />
        <Button
          class="p-button-text"
          icon="pi pi-arrow-right"
          label="next"
          @click="selectIndex(1)"
        />

        <Button class="ml-auto p-button-text" label="Cancel" @click="cancel" />
        <Button :disabled="!hasChanges" label="Save" @click="save" />
      </div>
    </template>
  </Dialog>
</template>
