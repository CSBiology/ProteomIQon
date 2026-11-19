function element(tagName, properties, children) {
	return { type: 'element', tagName, properties, children };
}

/** Keeps the Markdown source plain while labelling table fields for phones. */
export function responsiveTables() {
	return {
		name: 'responsive-tables-v1',
		element: {
			filter: ['table'],
			visit(table, ctx) {
				const thead = table.children.find((node) => node.type === 'element' && node.tagName === 'thead');
				const headerRow = thead?.children.find((node) => node.type === 'element' && node.tagName === 'tr');
				const labels = headerRow?.children.filter((node) => node.type === 'element' && node.tagName === 'th').map((node) => ctx.textContent(node).trim()) ?? [];
				if (!labels.length) return;

				function prepare(node) {
					if (node.type !== 'element') return node;
					const properties = { ...node.properties };
					const children = (node.children ?? []).map(prepare);
					if (node.tagName === 'thead' || node.tagName === 'tbody' || node.tagName === 'tfoot') properties.role = 'rowgroup';
					if (node.tagName === 'th') {
						properties.scope = 'col';
						properties.role = 'columnheader';
					}
					if (node.tagName === 'tr') {
						properties.role = 'row';
						let column = 0;
						for (const cell of children) {
							if (cell.type !== 'element' || !['td', 'th'].includes(cell.tagName)) continue;
							if (cell.tagName === 'td') cell.properties = { ...cell.properties, role: 'cell', 'data-label': labels[column] ?? '' };
							column += Number(cell.properties?.colSpan ?? 1);
						}
					}
					return { ...node, properties, children };
				}

				const prepared = prepare(table);
				const classes = prepared.properties.className ?? [];
				prepared.properties.className = [...(Array.isArray(classes) ? classes : String(classes).split(/\s+/)), 'responsive-table'];
				prepared.properties.role = 'table';
				const parentClasses = ctx.parent(table)?.properties?.className ?? [];
				if (parentClasses.includes('tool-io-table-scroll')) return prepared;
				return element('div', {
					className: ['content-table-scroll'],
					tabIndex: 0,
					role: 'region',
					ariaLabel: `${labels.join(', ')} table`,
				}, [prepared]);
			},
		},
	};
}
