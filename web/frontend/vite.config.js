import { defineConfig } from "vite";

// Dev server proxies /api to the FastAPI backend so the browser can call it
// without CORS. The backend runs on 127.0.0.1:8141 by default.
export default defineConfig({
  server: {
    port: 5173,
    proxy: {
      "/api": {
        target: process.env.MUSIC_API_TARGET || "http://127.0.0.1:8141",
        changeOrigin: true,
      },
    },
  },
  build: {
    outDir: "dist",
    emptyOutDir: true,
  },
});
