<script setup>
import { ref, shallowRef, watch } from 'vue'
import { useSessionStore } from '../store/session'

const session = useSessionStore()
const ketcher = shallowRef(null)
const visible = ref(false)

watch(
  () => session.sketchSmiles,
  value => {
    if (value === null) return
    visible.value = true
  }
)

function load() {
  const iframe = document.getElementById('ketcher')
  iframe.onload = async () => {
    while (!iframe.contentWindow.ketcher) {
      await new Promise(resolve => setTimeout(resolve, 100))
    }
    ketcher.value = iframe.contentWindow.ketcher
    ketcher.value.setMolecule(session.sketchSmiles)
  }
}

async function save() {
  const smiles = await ketcher.value.getSmiles()
  session.addSmiles(smiles)
  visible.value = false
}
</script>

<template>
  <Button
    class="mx-2 p-button-outlined"
    icon="pi pi-pencil"
    label="New Sketch"
    @click="session.newSketch()"
  />

  <Dialog
    v-model:visible="visible"
    modal
    @show="load"
    @hide="session.clearSketch()"
  >
    <template #header>Sketch Molecule</template>
    <iframe id="ketcher" src="/ketcher/index.html" />

    <template #footer>
      <div class="flex align-items-bottom text-left">
        <span class="p-float-label">
          <Dropdown
            v-model="session.sketchAlignTo"
            :options="session.complexes"
            option-label="name"
            option-value="index"
            class="w-10rem p-inputwrapper-filled"
            placeholder="none"
            show-clear
          />
          <label>Align to</label>
        </span>

        <label class="mx-4 flex align-items-center">
          <i class="mr-2 pi pi-plus-circle" />
          <span class="mr-2">Add Hydrogens</span>
          <InputSwitch v-model="session.sketchHydrogens" />
        </label>

        <div class="ml-auto">
          <Button
            class="p-button-text"
            icon="pi pi-times"
            label="Cancel"
            @click="visible = false"
          />
          <Button icon="pi pi-save" label="Save" @click="save" />
        </div>
      </div>
    </template>
  </Dialog>
</template>

<style>
#ketcher {
  border: 0;
  width: 90vw;
  height: calc(90vh - 12rem);
}
</style>
