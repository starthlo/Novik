import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';
import tailwindcss from '@tailwindcss/vite';
import dts from 'vite-plugin-dts';
import svgr from 'vite-plugin-svgr';

// https://vite.dev/config/
export default defineConfig({
  plugins: [react(), tailwindcss(), svgr(), dts()],
  // base: './',
  server: {
    port: 12000,
    strictPort: true,
    host: '0.0.0.0',
    proxy: {
      '/api': 'http://127.0.0.1:8000',
    },
    allowedHosts: ['novik.ai', '127.0.0.1', 'dev1.infractura.com'],
  },
  build: {
    target: ['es2015'],
    rollupOptions: {
      output: {
        manualChunks(id) {
          if (id.includes('react')) return 'react';
          if (id.includes('react-dom')) return 'react-dom';
          if (id.includes('react-router-dom')) return 'react-router';
          if (id.includes('@mui') || id.includes('tailwind')) return 'ui';
          if (id.includes('node_modules')) return 'vendor';
        },
      },
    },
  },
});
