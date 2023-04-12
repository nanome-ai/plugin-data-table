<script setup>
import { ref, shallowRef, watch } from 'vue'
import { useSessionStore } from '../store/session'

const session = useSessionStore()
const ketcher = shallowRef(null)

const visible = ref(false)
const align_to = ref(null)
const hydrogens = ref(false)

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
  session.addSmiles(smiles, align_to.value, hydrogens.value)
  visible.value = false
}
</script>

<template>
  <Button
    class="mx-2 p-button-outlined"
    icon="pi pi-pencil"
    label="New Sketch"
    @click="session.sketchSmiles = ''"
  />

  <Dialog
    v-model:visible="visible"
    modal
    @show="load"
    @hide="session.sketchSmiles = null"
  >
    <template #header>Sketch Molecule</template>
    <iframe id="ketcher" src="/ketcher/index.html" />

    <template #footer>
      <div class="flex align-items-bottom">
        <Dropdown
          v-model="align_to"
          :options="session.complexes"
          option-label="name"
          option-value="index"
          class="w-10rem p-inputwrapper-filled"
          placeholder="align to"
          show-clear
        />

        <label class="mx-4 flex align-items-center">
          <i class="mr-2 pi pi-plus-circle" />
          <span class="mr-2">Add Hydrogens</span>
          <InputSwitch v-model="hydrogens" />
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
  height: 70vh;
}
</style>
