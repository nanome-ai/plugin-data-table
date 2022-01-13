import { defineConfig } from 'vite'
import vue from '@vitejs/plugin-vue'

// https://vitejs.dev/config/
export default defineConfig({
  plugins: [vue()],
  server: {
    host: true,
    proxy: {
      '/ws': {
        target: 'http://localhost',
        changeOrigin: true,
        secure: false,
        ws: true
      }
    }
  }
})
