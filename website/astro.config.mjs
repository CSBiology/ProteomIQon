import { fileURLToPath } from 'node:url';
import { defineConfig } from 'astro/config';
import { satteri } from '@astrojs/markdown-satteri';
import react from '@astrojs/react';
import tailwindcss from '@tailwindcss/vite';
import { contentLinks } from './src/plugins/content-links.mjs';
import { toolContent } from './src/plugins/tool-content.mjs';
import { responsiveTables } from './src/plugins/responsive-tables.mjs';
import { alignmentReadability } from './src/plugins/alignment-readability.mjs';
import { mathTypesetting } from './src/plugins/math-typesetting.mjs';

const base = '/ProteomIQon';

export default defineConfig({
  site: 'https://csbiology.github.io',
  base,
  integrations: [react()],
  markdown: {
    processor: satteri({
      features: { smartPunctuation: false, math: true },
      hastPlugins: [contentLinks({ base, contentDir: fileURLToPath(new URL('./src/content', import.meta.url)) }), toolContent(), alignmentReadability(), responsiveTables(), mathTypesetting()],
    }),
  },
  vite: { plugins: [tailwindcss()] },
});
