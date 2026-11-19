import katex from 'katex';

/** Typesets Markdown math at build time, including an accessible MathML version. */
export function mathTypesetting() {
	return {
		name: 'math-typesetting-v1',
		element: {
			filter: ['pre', 'code'],
			visit(node, ctx) {
				const code = node.tagName === 'pre'
					? node.children.find((child) => child.type === 'element' && child.tagName === 'code')
					: node;
				const classes = code?.properties?.className ?? [];
				const display = node.tagName === 'pre' && classes.includes('math-display');
				if (!code || (!display && !classes.includes('math-inline'))) return;

				const rendered = {
					type: 'raw',
					value: katex.renderToString(ctx.textContent(code), {
						displayMode: display,
						output: 'htmlAndMathml',
						throwOnError: false,
					}),
				};
				if (!display) return rendered;
				return {
					type: 'element',
					tagName: 'div',
					properties: { className: ['tool-math-block'], tabIndex: 0, role: 'region', ariaLabel: 'Mathematical formula' },
					children: [rendered],
				};
			},
		},
	};
}
