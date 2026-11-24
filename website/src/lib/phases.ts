// The processing stages of the ProteomIQon pipeline, in pipeline order.
// Tool frontmatter refers to a stage by its id. Number, label and colour all come from here.
export const phases = [
	{ id: 'preprocessing', number: '01', label: 'Preprocessing', detail: 'mzML · mzlite · peptide database', tone: 'sky' },
	{ id: 'peptide-identification', number: '02', label: 'Peptide Identification', detail: 'spectrum matching · FDR', tone: 'orange' },
	{ id: 'peptide-quantification', number: '03', label: 'Peptide Quantification', detail: 'XIC · intensity', tone: 'teal' },
	{ id: 'protein-identification', number: '04', label: 'Protein Identification', detail: 'protein groups', tone: 'red' },
	{ id: 'alignment', number: '05', label: 'Alignment', detail: 'run alignment · peptide transfer', tone: 'violet' },
	{ id: 'protein-quantification', number: '06', label: 'Protein Quantification', detail: 'protein intensity · ratios', tone: 'green' },
] as const;

export type Phase = (typeof phases)[number];
export type PhaseId = Phase['id'];
export type PhaseTone = Phase['tone'];

export const phaseIds = phases.map((phase) => phase.id) as [PhaseId, ...PhaseId[]];

export function getPhase(id: PhaseId): Phase {
	return phases.find((phase) => phase.id === id)!;
}
