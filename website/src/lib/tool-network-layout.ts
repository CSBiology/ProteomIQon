type WorkflowNode = {
	id: number;
	tool: string;
	position: readonly [number, number];
	variant?: string;
	mobilityOnly?: boolean;
};

// Positions are percentages of the graph canvas. Edges are defined independently.
export const workflowNodes: readonly WorkflowNode[] = [
	{ id: 1, tool: 'peptide-db', position: [6, 28] },
	{ id: 2, tool: 'mzml-to-mzlite', position: [6, 57], variant: 'mzml-to-mzlite-ion-mobility' },
	{ id: 3, tool: 'peptide-spectrum-matching', position: [17, 43], variant: 'peptide-spectrum-matching-tims' },
	{ id: 4, tool: 'ms-fragger-to-psm', position: [27, 25], mobilityOnly: true },
	{ id: 5, tool: 'psm-statistics', position: [29, 59] },
	{ id: 6, tool: 'psm-based-quantification', position: [38, 39], variant: 'psm-based-quantification-tims' },
	{ id: 7, tool: 'protein-inference', position: [49, 61] },
	{ id: 8, tool: 'quant-based-alignment', position: [55.5, 19] },
	{ id: 9, tool: 'alignment-based-quantification', position: [64, 33] },
	{ id: 10, tool: 'alignment-based-quant-statistics', position: [73, 46] },
	{ id: 11, tool: 'add-deduced-peptides', position: [81, 24] },
	{ id: 12, tool: 'join-quant-with-peptide-ions', position: [81.4, 69.5] },
	{ id: 13, tool: 'label-free-protein-quantification', position: [88.1, 87.5] },
	{ id: 14, tool: 'labeled-protein-quantification', position: [88.1, 51.5] },
	{ id: 15, tool: 'ratio-lfq', position: [94.8, 69.5] },
];

export const workflowEdges = [
	[1, 3], [2, 3], [4, 6], [3, 5], [5, 6], [5, 7], [6, 8], [7, 11], [7, 12],
	[8, 9], [9, 10], [10, 11], [11, 12], [12, 13], [12, 14], [14, 15],
] as const;

type Circle = { x: number; y: number; radius: number };

export function connectToolCircles(from: Circle, to: Circle) {
	const dx = to.x - from.x;
	const dy = to.y - from.y;
	const distance = Math.hypot(dx, dy);
	const ux = dx / distance;
	const uy = dy / distance;
	return `M ${from.x + ux * from.radius} ${from.y + uy * from.radius} L ${to.x - ux * to.radius} ${to.y - uy * to.radius}`;
}
