import { createHash } from 'node:crypto';
import { readdirSync } from 'node:fs';
import { basename, join } from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';

const SKIP = /^([a-z][a-z0-9+.-]*:|\/|#|\?)/i;
const MD_LINK = /^([^?#]+\.md)([?#].*)?$/i;
const DIR_LINK = /^([^?#]*\/)([?#].*)?$/;
const KEBAB = /^[a-z0-9]+(-[a-z0-9]+)*$/;
// Bump after changing the plugin code.
const VERSION = '3';

/**
 * Rewrites relative links between Markdown files (`./peptide-db.md`, `../guides/help.md#x`)
 * to the routes of the pages built from them, including the site base.
 *
 * Links to the tools directory (`../tools/`) are rewritten as well. A relative `.md` link that
 * points to no known page is left as it is and logged as a `[content-links]` warning.
 *
 * The plugin name carries a hash of VERSION and the route map because the content layer cache
 * does not fingerprint plugin code or the set of files. Bump VERSION after changing the plugin,
 * and adding, renaming or removing a content file changes the hash by itself.
 */
export function contentLinks({ base, contentDir }) {
	const routes = new Map();
	const register = (dir, prefix) => {
		for (const file of readdirSync(join(contentDir, dir))) {
			if (!file.endsWith('.md')) continue;
			const name = basename(file, '.md');
			if (!KEBAB.test(name)) {
				throw new Error(`Content file "${dir}/${file}" must have a lowercase kebab-case name (for example "peptide-db.md"), because file names become URLs.`);
			}
			routes.set(pathToFileURL(join(contentDir, dir, file)).href, `${prefix}/${name}/`);
		}
	};
	register('tools', '/tools');
	register('guides', '/guide');
	const prefix = base.replace(/\/+$/, '');
	const toolsDir = pathToFileURL(join(contentDir, 'tools')).href;
	const dirRoutes = new Map([[toolsDir + '/', '/tools/'], [toolsDir, '/tools/']]);
	const hash = createHash('sha1')
		.update(VERSION + JSON.stringify([...routes].sort(([a], [b]) => (a < b ? -1 : a > b ? 1 : 0))))
		.digest('hex')
		.slice(0, 8);
	const lowerKeys = new Map([...routes.keys(), ...dirRoutes.keys()].map((key) => [key.toLowerCase(), key]));

	return {
		name: `content-links-${hash}`,
		element: {
			filter: ['a'],
			visit(node, ctx) {
				const href = node.properties?.href;
				if (typeof href !== 'string' || SKIP.test(href) || !ctx.fileURL) return;
				const md = href.match(MD_LINK);
				const dir = md ? null : href.match(DIR_LINK);
				const match = md ?? dir;
				if (!match) return;
				const target = new URL(match[1], ctx.fileURL).href;
				const route = md ? routes.get(target) : dirRoutes.get(target);
				if (!route) {
					const near = lowerKeys.get(target.toLowerCase());
					// Warn instead of throwing: the glob loader catches errors from plugins, drops the
					// whole page and still lets the build succeed, which is worse than one dead link.
					if (near) console.warn(`[content-links] Link "${href}" in ${fileURLToPath(ctx.fileURL)} differs only in letter case from an existing page. Fix the case.`);
					else if (md) console.warn(`[content-links] Link "${href}" in ${fileURLToPath(ctx.fileURL)} points to a Markdown file that is not a page.`);
					return;
				}
				ctx.setProperty(node, 'href', prefix + route + (match[2] ?? ''));
			},
		},
	};
}
