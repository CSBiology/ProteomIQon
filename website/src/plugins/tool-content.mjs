const SECTION_DETAILS = {
	inputs: {
		kind: 'inputs',
		icon: 'IN',
		summary: 'Files and command-line flags',
	},
	outputs: {
		kind: 'outputs',
		icon: 'OUT',
		summary: 'Generated data and log files',
	},
	'parameter file': {
		kind: 'parameters',
		icon: 'CFG',
		summary: 'Default values and tool settings',
	},
	parameters: {
		kind: 'parameters',
		icon: 'CFG',
		summary: 'Default values and tool settings',
	},
	'input file handling': {
		kind: 'inputs',
		icon: 'IN',
		summary: 'File selection and pairing',
	},
	'parameter generation': {
		kind: 'parameters',
		icon: 'CFG',
		summary: 'JSON and parameter-generation examples',
	},
};

function element(tagName, properties = {}, children = []) {
	return { type: 'element', tagName, properties, children };
}

function text(value) {
	return { type: 'text', value };
}

// Keep Markdown tables semantic and give each one its own scroll container on
// small screens, including the wavelet table nested inside a Markdown list.
function wrapTables(nodes, title) {
	return nodes.map((node) => {
		if (node.type !== 'element') return node;
		const prepared = { ...node, children: wrapTables(node.children ?? [], title) };
		if (node.tagName !== 'table') return prepared;
		return element('div', {
			className: ['tool-io-table-scroll'],
			tabIndex: 0,
			role: 'region',
			ariaLabel: `${title} table`,
		}, [prepared]);
	});
}

function accordion(heading, title, nodes) {
	const details = SECTION_DETAILS[title.toLowerCase()];
	const id = heading.properties?.id ?? title.toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/^-|-$/g, '');
	return element('details', { id, className: ['accordion', 'tool-io-accordion', `tool-io-${details.kind}`] }, [
		element('summary', { className: ['accordion-header'] }, [
			element('span', { className: ['accordion-icon'], ariaHidden: 'true' }, [text(details.icon)]),
			element('span', { className: ['accordion-heading'] }, [
				element('strong', { role: 'heading', ariaLevel: 2 }, [text(title)]),
				element('small', {}, [text(details.summary)]),
			]),
			element('span', { className: ['accordion-chevron'], ariaHidden: 'true' }),
		]),
		element('div', { className: ['accordion-content'] }, wrapTables(nodes, title)),
	]);
}

function addClass(node, className) {
	const classes = node.properties?.className ?? [];
	return { ...node, properties: { ...node.properties, className: [...(Array.isArray(classes) ? classes : String(classes).split(/\s+/)), className] } };
}

function workflowList(list) {
	const children = list.children.map((child) => {
		if (child.type !== 'element' || child.tagName !== 'li') return child;
		const content = child.children.map((node) => {
			if (node.type !== 'element' || node.tagName !== 'p') return node;
			const meaningful = node.children.filter((part) => part.type !== 'text' || part.value.trim());
			return meaningful.length === 1 && meaningful[0].type === 'element' && meaningful[0].tagName === 'strong'
				? addClass(node, 'tool-step-heading') : node;
		});
		return { ...child, children: [element('div', { className: ['tool-step-copy'] }, content)] };
	});
	const prepared = addClass(list, 'tool-explanation-steps');
	return { ...prepared, properties: { ...prepared.properties, role: 'list' }, children };
}

/** Applies shared presentation to tool Markdown without adding markup to content files. */
export function toolContent() {
	return {
		name: 'tool-content-v3',
		after(root, ctx) {
			if (!/\/content\/tools\/[^/]+\.md$/.test(String(ctx.fileURL ?? ''))) return;

			const children = root.children ?? [];
			const nextChildren = [];

			for (let i = 0; i < children.length; i += 1) {
				const child = children[i];
				const title = child.type === 'element' && child.tagName === 'h2' ? ctx.textContent(child).trim() : '';

				if (!title) {
					nextChildren.push(child);
					continue;
				}

				let end = i + 1;
				while (end < children.length && !(children[end].type === 'element' && children[end].tagName === 'h2')) {
					end += 1;
				}

				const content = children.slice(i + 1, end);
				if (Object.hasOwn(SECTION_DETAILS, title.toLowerCase())) {
					nextChildren.push(accordion(child, title, content));
				} else if (/^how does\b|^what['’]s happening|^model, calibration and q-values$/i.test(title)) {
					nextChildren.push(addClass(child, 'tool-workflow-heading'));
					nextChildren.push(...content.map((node) => node.type === 'element' && node.tagName === 'ol' ? workflowList(node) : node));
				} else {
					nextChildren.push(child, ...content);
				}
				i = end - 1;
			}

			ctx.replaceNode(root, { type: 'root', children: nextChildren });
		},
	};
}
