import { getCollection } from 'astro:content';
import { phases } from './phases';
import { withBase } from './site';

const phaseIndex = (id: string) => phases.findIndex((phase) => phase.id === id);

export async function getSortedTools() {
	const tools = await getCollection('tools');
	return tools.sort((a, b) => phaseIndex(a.data.phase) - phaseIndex(b.data.phase) || a.data.order - b.data.order || a.id.localeCompare(b.id));
}

export const toolHref = (id: string) => withBase(`/tools/${id}/`);
