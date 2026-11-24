type Point = { x: number; y: number };
type Circle = Point & { radius: number };
type Size = { width: number; height: number };
type Box = Point & Size;
type HalfPlane = { normal: Point; limit: number };
type RegionGroup = { circles: Circle[]; label: Size };

const dot = (point: Point, normal: Point) => point.x * normal.x + point.y * normal.y;

export function layoutCategoryRegions(canvas: Size, groups: RegionGroup[], edges: [Point, Point][]) {
	const gutter = 5;
	const cells: HalfPlane[][] = groups.map(() => [
		{ normal: { x: -1, y: 0 }, limit: -gutter },
		{ normal: { x: 1, y: 0 }, limit: canvas.width - gutter },
		{ normal: { x: 0, y: -1 }, limit: -gutter },
		{ normal: { x: 0, y: 1 }, limit: canvas.height - gutter },
	]);
	for (const [i, group] of groups.entries()) {
		for (let j = i + 1; j < groups.length; j++) {
			const { normal, middle, gap } = separatingPlane(group.circles, groups[j].circles);
			const clearance = Math.min(6, gap / 4);
			cells[i].push({ normal, limit: middle - clearance });
			cells[j].push({ normal: { x: -normal.x, y: -normal.y }, limit: -middle - clearance });
		}
	}

	const allCircles = groups.flatMap(({ circles }) => circles);
	return groups.map(({ circles, label }, index) => {
		const cell = cells[index];
		const labelPosition = findLabel(canvas, circles, allCircles, edges, label, cell);
		const outline = circles.flatMap((circle) => Array.from({ length: 24 }, (_, i) => {
			const angle = Math.PI * i / 12;
			const padding = Math.min(38, circle.radius * .62);
			const radius = circle.radius + padding + 5 * Math.sin(3 * angle + index) + 3 * Math.cos(2 * angle - index);
			return { x: circle.x + Math.cos(angle) * radius, y: circle.y + Math.sin(angle) * radius };
		}));
		if (labelPosition) outline.push(...corners({ ...label, ...labelPosition }, 12));

		const polygon = cell.reduce(clipPolygon, convexHull(outline));
		return { path: roundPolygon(polygon), labelPosition };
	});
}

// Prefer upright separators along the horizontal workflow, turning only as needed to clear the tools.
function separatingPlane(first: Circle[], second: Circle[]) {
	const planes = Array.from({ length: 90 }, (_, i) => {
		const angle = i * Math.PI / 45;
		const normal = { x: Math.cos(angle), y: Math.sin(angle) };
		const a = Math.max(...first.map((circle) => dot(circle, normal) + circle.radius));
		const b = Math.min(...second.map((circle) => dot(circle, normal) - circle.radius));
		return { normal, middle: (a + b) / 2, gap: b - a };
	});
	const widest = planes.reduce((a, b) => a.gap > b.gap ? a : b);
	const minimumGap = Math.max(12, widest.gap * .65);
	return planes.reduce((best, plane) => (
		plane.gap >= minimumGap && plane.normal.x > best.normal.x ? plane : best
	), widest);
}

function corners({ x, y, width, height }: Box, padding: number) {
	return [
		{ x: x - padding, y: y - padding },
		{ x: x + width + padding, y: y - padding },
		{ x: x + width + padding, y: y + height + padding },
		{ x: x - padding, y: y + height + padding },
	];
}

function findLabel(canvas: Size, group: Circle[], circles: Circle[], edges: [Point, Point][], label: Size, cell: HalfPlane[]) {
	const left = Math.min(...group.map((circle) => circle.x - circle.radius));
	const right = Math.max(...group.map((circle) => circle.x + circle.radius));
	const top = Math.min(...group.map((circle) => circle.y - circle.radius));
	const bottom = Math.max(...group.map((circle) => circle.y + circle.radius));
	const preferred = { x: left + 8, y: top - label.height - 12 };
	let best: Point | undefined;
	let score = Infinity;

	for (let y = Math.max(18, top - 84); y <= Math.min(canvas.height - label.height - 18, bottom + 46); y += 6) {
		for (let x = Math.max(18, left - 60); x <= Math.min(canvas.width - label.width - 18, right + 48); x += 6) {
			const distance = Math.hypot(x - preferred.x, y - preferred.y);
			if (distance >= score) continue;

			const box = { x, y, ...label };
			if (corners(box, 12).some((point) => cell.some(({ normal, limit }) => dot(point, normal) > limit))) continue;
			if (circles.some((circle) => boxTouchesCircle(box, circle))) continue;
			if (edges.some(([from, to]) => lineTouchesBox(from, to, box))) continue;

			best = { x, y };
			score = distance;
		}
	}
	return best;
}

function boxTouchesCircle(box: Box, circle: Circle) {
	const x = Math.max(box.x, Math.min(circle.x, box.x + box.width));
	const y = Math.max(box.y, Math.min(circle.y, box.y + box.height));
	return Math.hypot(circle.x - x, circle.y - y) < circle.radius + 9;
}

function lineTouchesBox(from: Point, to: Point, box: Box) {
	let enter = 0;
	let leave = 1;
	for (const [start, end, min, max] of [
		[from.x, to.x, box.x - 6, box.x + box.width + 6],
		[from.y, to.y, box.y - 6, box.y + box.height + 6],
	]) {
		const delta = end - start;
		if (delta === 0) {
			if (start < min || start > max) return false;
			continue;
		}
		const a = (min - start) / delta;
		const b = (max - start) / delta;
		enter = Math.max(enter, Math.min(a, b));
		leave = Math.min(leave, Math.max(a, b));
	}
	return enter <= leave;
}

function convexHull(points: Point[]) {
	const sorted = [...points].sort((a, b) => a.x - b.x || a.y - b.y);
	const half = (ordered: Point[]) => {
		const hull: Point[] = [];
		for (const point of ordered) {
			while (hull.length >= 2) {
				const b = hull[hull.length - 1];
				const a = hull[hull.length - 2];
				if ((b.x - a.x) * (point.y - a.y) > (b.y - a.y) * (point.x - a.x)) break;
				hull.pop();
			}
			hull.push(point);
		}
		return hull.slice(0, -1);
	};
	return [...half(sorted), ...half(sorted.reverse())];
}

function clipPolygon(polygon: Point[], { normal, limit }: HalfPlane) {
	return polygon.flatMap((a, i) => {
		const b = polygon[(i + 1) % polygon.length];
		const da = dot(a, normal) - limit;
		const db = dot(b, normal) - limit;
		const points = da <= 0 ? [a] : [];
		if ((da <= 0) !== (db <= 0)) {
			const t = da / (da - db);
			points.push({ x: a.x + (b.x - a.x) * t, y: a.y + (b.y - a.y) * t });
		}
		return points;
	});
}

function roundPolygon(polygon: Point[]) {
	const toward = (from: Point, to: Point) => {
		const distance = Math.hypot(to.x - from.x, to.y - from.y);
		const t = distance ? Math.min(1 / 3, 24 / distance) : 0;
		return `${from.x + (to.x - from.x) * t} ${from.y + (to.y - from.y) * t}`;
	};
	return polygon.map((point, i) => {
		const before = polygon[(i + polygon.length - 1) % polygon.length];
		const after = polygon[(i + 1) % polygon.length];
		return `${i ? 'L' : 'M'} ${toward(point, before)} Q ${point.x} ${point.y} ${toward(point, after)}`;
	}).join(' ') + ' Z';
}
