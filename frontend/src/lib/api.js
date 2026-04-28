/**
 * primerdesignr API client
 */

const API_BASE = import.meta.env.VITE_API_URL || (import.meta.env.DEV ? 'http://localhost:8000' : '');

function apiUrl(path) {
	return `${API_BASE}${path}`;
}

function normalizeSequence(text) {
	return text
		.split('\n')
		.filter((line) => !line.trim().startsWith('>'))
		.join('')
		.replace(/\s+/g, '')
		.toUpperCase();
}

function csvCell(value) {
	const text = String(value ?? '');
	if (/[",\n\r]/.test(text)) return `"${text.replace(/"/g, '""')}"`;
	return text;
}

/**
 * @typedef {Object} Conditions
 * @property {number} na_mm
 * @property {number} mg_mm
 * @property {number} dna_nm
 */

/**
 * @typedef {Object} AnalyzeRequest
 * @property {Record<string, string>} primers
 * @property {Conditions} conditions
 * @property {boolean} include_cross_dimers
 */

/**
 * Parse raw text input into a primers dict.
 * Handles: "NAME SEQ", "NAME\tSEQ", ">NAME\nSEQ", bare "SEQ"
 * @param {string} raw
 * @returns {Record<string, string>}
 */
export function parsePrimerInput(raw) {
	/** @type {Record<string, string>} */
	const primers = {};
	const lines = raw.trim().split('\n');
	let counter = 1;

	for (let i = 0; i < lines.length; i++) {
		const line = lines[i].trim();
		if (!line) continue;

		// FASTA
		if (line.startsWith('>')) {
			const name = line.slice(1).trim();
			const seqLines = [];
			while (i + 1 < lines.length && !lines[i + 1].trim().startsWith('>')) {
				i++;
				if (lines[i].trim()) seqLines.push(lines[i].trim());
			}
			const seq = normalizeSequence(seqLines.join('\n'));
			if (name && /^[ATGC]+$/.test(seq) && seq.length >= 10) {
				primers[name] = seq;
			}
			continue;
		}

		// Tab, comma, or space separated
		const match = line.match(/^(.+?)[\t, ]+([ATGCatgc]{10,})$/);
		if (match) {
			primers[match[1].trim()] = match[2].toUpperCase();
			continue;
		}

		// Bare sequence
		const seq = line.toUpperCase();
		if (/^[ATGC]+$/.test(seq) && seq.length >= 10) {
			primers[`Primer_${counter++}`] = seq;
		}
	}

	return primers;
}

/**
 * Call the /analyze endpoint
 * @param {Record<string, string>} primers
 * @param {Conditions} conditions
 * @returns {Promise<Object>}
 */
export async function analyzePrimers(primers, conditions = {}) {
	const res = await fetch(apiUrl('/analyze'), {
		method: 'POST',
		headers: { 'Content-Type': 'application/json' },
		body: JSON.stringify({
			primers,
			conditions: {
				na_mm: conditions.na_mm ?? 50,
				mg_mm: conditions.mg_mm ?? 0,
				dntp_mm: conditions.dntp_mm ?? 0,
				dna_nm: conditions.dna_nm ?? 250,
			},
			include_cross_dimers: true,
		}),
	});

	if (!res.ok) {
		const err = await res.json().catch(() => ({ detail: res.statusText }));
		throw new Error(err.detail || 'Analysis failed');
	}

	return res.json();
}

/**
 * Design PCR primers from a DNA template.
 * @param {Object} payload
 * @returns {Promise<Object>}
 */
export async function designPrimers(payload) {
	const res = await fetch(apiUrl('/design'), {
		method: 'POST',
		headers: { 'Content-Type': 'application/json' },
		body: JSON.stringify(payload),
	});

	if (!res.ok) {
		const err = await res.json().catch(() => ({ detail: res.statusText }));
		throw new Error(err.detail || 'Primer design failed');
	}

	return res.json();
}

/**
 * Validate Golden Gate overhangs.
 * @param {string[]} overhangs
 * @param {string} enzyme
 * @param {Record<string, string>|null} sequences
 * @returns {Promise<Object>}
 */
export async function checkGoldenGate(overhangs, enzyme = 'BsaI', sequences = null) {
	const res = await fetch(apiUrl('/golden-gate'), {
		method: 'POST',
		headers: { 'Content-Type': 'application/json' },
		body: JSON.stringify({ overhangs, enzyme, sequences }),
	});

	if (!res.ok) {
		const err = await res.json().catch(() => ({ detail: res.statusText }));
		throw new Error(err.detail || 'Golden Gate check failed');
	}

	return res.json();
}

/**
 * Export results as CSV string
 * @param {Object} results - AnalyzeResponse from the API
 * @returns {string}
 */
export function exportCSV(results) {
	const rows = [
		['Name', 'Sequence', 'Length', 'GC%', 'Tm (°C)',
		 'Hairpin ΔG (SL)', 'Hairpin ΔG (MW)', 'Self-dimer ΔG',
		 'Dot-bracket', 'Warnings'].join(','),
	];

	for (const [name, p] of Object.entries(results.primers)) {
		rows.push([
			csvCell(name),
			csvCell(p.sequence),
			p.length,
			csvCell((p.gc_content * 100).toFixed(1) + '%'),
			p.tm,
			p.hairpin.dg_santalucia,
			p.hairpin.dg_mathews,
			p.homodimer.dg,
			csvCell(p.hairpin.dot_bracket),
			csvCell(p.warnings.join('; ')),
		].join(','));
	}

	if (results.cross_dimers?.length) {
		rows.push('');
		rows.push('Cross-dimer pairs');
		rows.push('Primer A,Primer B,ΔG (kcal/mol),Problematic');
		for (const cd of results.cross_dimers) {
			rows.push([csvCell(cd.primer_a), csvCell(cd.primer_b), cd.dg, cd.is_problematic].join(','));
		}
	}

	return rows.join('\n');
}

export function exportDesignCSV(results) {
	const rows = [
		[
			'Rank',
			'Forward',
			'Forward sequence',
			'Reverse',
			'Reverse sequence',
			'Product size',
			'Tm difference',
			'Heterodimer ΔG',
			'Warnings',
		].join(','),
	];

	for (const c of results.candidates ?? []) {
		rows.push([
			c.rank,
			'Forward',
			csvCell(c.forward.sequence),
			'Reverse',
			csvCell(c.reverse.sequence),
			c.product_size,
			c.tm_difference,
			c.heterodimer.dg,
			csvCell(c.warnings.join('; ')),
		].join(','));
	}

	return rows.join('\n');
}

export { normalizeSequence };
