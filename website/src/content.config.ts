import { defineCollection, reference } from 'astro:content';
import { glob } from 'astro/loaders';
import { z } from 'astro/zod';
import { phaseIds } from './lib/phases';

const tools = defineCollection({
	loader: glob({ pattern: '*.md', base: './src/content/tools' }),
	schema: z.object({
		name: z.string(),
		phase: z.enum(phaseIds),
		subgroup: z.string().min(1).optional(),
		order: z.number().int().positive(),
		summary: z.string().min(10),
		inputsFrom: z.array(reference('tools')).default([]),
	}),
});

const guides = defineCollection({
	loader: glob({ pattern: '*.md', base: './src/content/guides' }),
	schema: ({ image }) =>
		z.object({
			title: z.string(),
			eyebrow: z.string(),
			intro: z.string(),
			description: z.string(),
			order: z.number().int(),
			contributors: z.array(z.object({ name: z.string(), role: z.string().optional(), url: z.url().optional(), image: image().optional() })).optional(),
		}),
});

export const collections = { tools, guides };
