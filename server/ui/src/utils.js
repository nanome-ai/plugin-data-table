export const randStr = n =>
  Array.from({ length: n }, () =>
    Math.floor(36 * Math.random()).toString(36)
  ).join('')
