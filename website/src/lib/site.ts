const configuredBase = import.meta.env.BASE_URL;

export const baseUrl = configuredBase.endsWith('/') ? configuredBase : `${configuredBase}/`;

export function withBase(path: string) {
	if (path === '/') return baseUrl;
	if (path.startsWith(baseUrl)) return path;
	return `${baseUrl}${path.replace(/^\/+/, '')}`;
}
