import { computed } from 'vue'

export function useMux(conditionFn, trueRef, falseRef) {
  return computed({
    get: () => (conditionFn() ? trueRef.value : falseRef.value),
    set: value => {
      const target = conditionFn() ? trueRef : falseRef
      target.value = value
    }
  })
}

export function useVModel(props, emit, name = 'modelValue') {
  return computed({
    get: () => props[name],
    set: value => emit(`update:${name}`, value)
  })
}
