<script setup>
import { computed, reactive, ref } from 'vue'
import { useSessionStore } from '../store/session'

const operators = [
  { name: '+ add', fn: (a, b) => a + b, arity: 2 },
  { name: '- subtract', fn: (a, b) => a - b, arity: 2 },
  { name: '× multiply', fn: (a, b) => a * b, arity: 2 },
  { name: '÷ divide', fn: (a, b) => a / b, arity: 2 },
  { name: '^ pow', fn: (a, b) => Math.pow(a, b), arity: 2 },
  { name: 'ln', fn: a => Math.log(a), arity: 1 },
  { name: 'log2', fn: a => Math.log2(a), arity: 1 },
  { name: 'log10', fn: a => Math.log10(a), arity: 1 },
  { name: 'sin', fn: a => Math.sin(a), arity: 1 },
  { name: 'cos', fn: a => Math.cos(a), arity: 1 },
  { name: 'tan', fn: a => Math.tan(a), arity: 1 }
]

const session = useSessionStore()
const menu = ref(null)

const data = reactive({
  name: '',
  operator: null,
  isOperand1Column: true,
  isOperand2Column: false,
  column1: null,
  column2: null,
  value1: null,
  value2: null
})

const canAdd = computed(() => {
  if (!data.name || !data.operator) return false
  const aValue = data.isOperand1Column ? data.column1 : data.value1
  const bValue = data.isOperand2Column ? data.column2 : data.value2
  if (aValue === null) return false
  if (data.operator.arity === 2 && bValue === null) return false
  return true
})

const addColumn = () => {
  const values = Array.from({ length: session.frames.length }, () => '')
  for (let i = 0; i < session.frames.length; i++) {
    const f = session.frames[i]
    const a = data.isOperand1Column ? f[data.column1] : data.value1
    const b = data.isOperand2Column ? f[data.column2] : data.value2
    values[i] = data.operator.fn(+a, +b)
    if (isNaN(values[i])) {
      values[i] = ''
    } else {
      values[i] = Math.round(values[i] * 1000) / 1000
    }
  }

  session.addColumn(data.name, values)
  menu.value.toggle(false)

  data.name = ''
  data.operator = null
  data.column1 = null
  data.column2 = null
  data.value1 = null
  data.value2 = null
}
</script>

<template>
  <Button
    class="p-button-outlined"
    label="new column"
    icon="pi pi-plus"
    @click="e => menu.toggle(e)"
  />

  <OverlayPanel ref="menu">
    <div class="mt-3 w-15rem flex flex-column gap-5">
      <span class="p-float-label">
        <InputText
          v-model="data.name"
          class="w-full p-inputwrapper-filled"
          placeholder="enter a name"
        />
        <label>Name</label>
      </span>

      <span class="p-float-label">
        <Dropdown
          v-model="data.operator"
          :options="operators"
          class="w-full p-inputwrapper-filled"
          option-label="name"
          placeholder="select an operator"
        />
        <label>Operator</label>
      </span>

      <template v-if="data.operator">
        <div>
          <span class="mt-3 mb-2 p-float-label">
            <Dropdown
              v-if="data.isOperand1Column"
              v-model="data.column1"
              :options="session.numericColumns"
              class="w-full p-inputwrapper-filled"
              placeholder="select a column"
            />
            <InputNumber
              v-else
              v-model="data.value1"
              class="w-full p-inputwrapper-filled"
              placeholder="enter a value"
            />
            <label>Operand 1</label>
          </span>

          <div class="field-checkbox text-sm">
            <Checkbox v-model="data.isOperand1Column" id="op1" binary />
            <label for="op1">Use column for operand 1</label>
          </div>
        </div>

        <div v-if="data.operator.arity === 2">
          <span class="mb-2 p-float-label">
            <Dropdown
              v-if="data.isOperand2Column"
              v-model="data.column2"
              :options="session.numericColumns"
              class="w-full p-inputwrapper-filled"
              placeholder="select a column"
            />
            <InputNumber
              v-else
              v-model="data.value2"
              class="w-full p-inputwrapper-filled"
              placeholder="enter a value"
            />
            <label>Operand 2</label>
          </span>

          <div class="field-checkbox text-sm">
            <Checkbox v-model="data.isOperand2Column" id="op2" binary />
            <label for="op2">Use column for operand 2</label>
          </div>
        </div>
      </template>

      <Button
        :disabled="!canAdd"
        label="create column"
        icon="pi pi-plus"
        class="p-button-outlined w-full"
        @click="addColumn"
      />
    </div>
  </OverlayPanel>
</template>
