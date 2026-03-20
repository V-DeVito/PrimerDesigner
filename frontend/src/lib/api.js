/**
 * primerdesignr API client
 */

const API_BASE = import.meta.env.VITE_API_URL || 'http://localhost:8000';

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
			i++;
			if (i < lines.length) {
				const seq = lines[i].trim().toUpperCase();
				if (/^[ATGC]+$/i.test(seq)) {
					primers[name] = seq;
				}
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
 * @param {Object} conditions
 * @returns {Promise<Object>}
 */
export async function analyzePrimers(primers, conditions = {}) {
	const res = await fetch(`${API_BASE}/analyze`, {
		method: 'POST',
		headers: { 'Content-Type': 'application/json' },
		body: JSON.stringify({
			primers,
			conditions: {
				na_mm: conditions.na_mm ?? 50,
				mg_mm: conditions.mg_mm ?? 0,
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
 * Call the /golden-gate endpoint
 * @param {string[]} overhangs
 * @param {string} enzyme
 * @param {Record<string, string>} [sequences]
 * @param {Record<string, string>} [primers]
 * @returns {Promise<Object>}
 */
export async function validateGoldenGate(overhangs, enzyme = 'BsaI', sequences, primers) {
	const body = { overhangs, enzyme };
	if (sequences) body.sequences = sequences;
	if (primers) body.primers = primers;

	const res = await fetch(`${API_BASE}/golden-gate`, {
		method: 'POST',
		headers: { 'Content-Type': 'application/json' },
		body: JSON.stringify(body),
	});

	if (!res.ok) {
		const err = await res.json().catch(() => ({ detail: res.statusText }));
		throw new Error(err.detail || 'Validation failed');
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
		 'Self-dimer 3\' ΔG', 'Dot-bracket', 'Warnings'].join(','),
	];

	for (const [name, p] of Object.entries(results.primers)) {
		rows.push([
			name,
			p.sequence,
			p.length,
			(p.gc_content * 100).toFixed(1) + '%',
			p.tm,
			p.hairpin.dg_santalucia,
			p.hairpin.dg_mathews,
			p.homodimer.dg,
			p.homodimer.dg_3prime ?? '',
			`"${p.hairpin.dot_bracket}"`,
			`"${p.warnings.join('; ')}"`,
		].join(','));
	}

	if (results.cross_dimers?.length) {
		rows.push('');
		rows.push('Cross-dimer pairs');
		rows.push('Primer A,Primer B,ΔG (kcal/mol),3\' ΔG,Problematic');
		for (const cd of results.cross_dimers) {
			rows.push(`${cd.primer_a},${cd.primer_b},${cd.dg},${cd.dg_3prime ?? ''},${cd.is_problematic}`);
		}
	}

	return rows.join('\n');
}
