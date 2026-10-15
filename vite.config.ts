import { defineConfig } from 'vite';
import { svelte } from '@sveltejs/vite-plugin-svelte';

export default defineConfig({
  base: './',
  plugins: [svelte()],
  server: {
    proxy: {
      '/preparic': { target: 'http://127.0.0.1:5813', ws: true },
      '/seamless-client.js': 'http://127.0.0.1:5813',
    },
  },
});
