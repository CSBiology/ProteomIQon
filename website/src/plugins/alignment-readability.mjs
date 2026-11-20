const ALIGNMENT_TOOLS = /\/tools\/(quant-based-alignment|alignment-based-quantification|alignment-based-quant-statistics)\.md$/;

function element(tagName, properties, children) {
	return { type: 'element', tagName, properties, children };
}

function addClass(node, name) {
	const classes = node.properties?.className ?? [];
	return { ...node, properties: { ...node.properties, className: [...(Array.isArray(classes) ? classes : String(classes).split(/\s+/)), name] } };
}

/** Adds figures, captions and reading cues without presentation markup in Markdown. */
export function alignmentReadability() {
	return {
		name: 'alignment-readability-v2',
		after(root, ctx) {
			if (!ALIGNMENT_TOOLS.test(String(ctx.fileURL ?? ''))) return;

			function decorate(nodes) {
				const result = [];
				for (let i = 0; i < nodes.length; i += 1) {
					const node = nodes[i];
					if (node.type !== 'element') { result.push(node); continue; }
					const meaningful = (node.children ?? []).filter((child) => child.type !== 'text' || child.value.trim());
					if (node.tagName === 'p' && meaningful.length === 1 && meaningful[0].type === 'element' && meaningful[0].tagName === 'img') {
						const image = meaningful[0];
						const overview = String(image.properties?.src ?? '').endsWith('/quantbasedalignment.png');
						const figure = element('figure', { className: ['alignment-figure', ...(overview ? ['alignment-figure-overview'] : [])] }, [image]);
						let captionIndex = i + 1;
						while (nodes[captionIndex]?.type === 'text' && !nodes[captionIndex].value.trim()) captionIndex += 1;
						const following = nodes[captionIndex];
						if (following?.type === 'element' && following.tagName === 'p' && /^(The figure shows|In the figure)\b/.test(ctx.textContent(following).trimStart())) {
							figure.children.push(element('figcaption', {}, following.children));
							i = captionIndex;
						}
						result.push(figure);
						continue;
					}
					let prepared = { ...node, children: decorate(node.children ?? []) };
					if (node.tagName === 'p' && /^(A peak is always chosen|The highest probability|A run with fewer than 200 training ions)/.test(ctx.textContent(node))) prepared = addClass(prepared, 'alignment-note');
					if (node.tagName === 'ul' && meaningful.length === 2 && meaningful.every((child) => child.type === 'element' && child.tagName === 'li') && /^AlignmentScore:/.test(ctx.textContent(meaningful[0])) && /^AlignmentQValue:/.test(ctx.textContent(meaningful[1]))) prepared = addClass(prepared, 'alignment-score-definitions');
					result.push(prepared);
				}
				return result;
			}
			ctx.replaceNode(root, { type: 'root', children: decorate(root.children ?? []) });
		},
	};
}
