<script>
	import { parsePrimerInput, analyzePrimers, exportCSV } from '$lib/api.js';

	// ── State (Svelte 5 runes) ─────────────────────────────
	let input = $state(`F_Chr\tGTCTTCACATCGGTTTGAAAGGAGG
R_Chr\tAACCCGCTCCGATTAAAGCTACTTT
F_TMM\tCCACAATGGAAGCAGATGTGATCGA
R_TMM\tACCACCCCAGTTATCCTGCTTTTCA
F_Kan\tTTTATTGATCTTGGGAGAAGCGGCA
R_Kan\tGTCCTGGGTTTCAAGCATTAGTCCA`);

	let naConc = $state(50);
	let mgConc = $state(0);
	let dnaConc = $state(250);
	let loading = $state(false);
	let error = $state('');
	let results = $state(null);

	// ── Derived ────────────────────────────────────────────
	let primerCount = $derived(Object.keys(parsePrimerInput(input)).length);
	let primerNames = $derived(results ? Object.keys(results.primers) : []);

	// ── Actions ────────────────────────────────────────────
	async function analyze() {
		const primers = parsePrimerInput(input);
		if (Object.keys(primers).length === 0) {
			error = 'No valid primer sequences found';
			return;
		}

		loading = true;
		error = '';
		results = null;

		try {
			results = await analyzePrimers(primers, {
				na_mm: naConc,
				mg_mm: mgConc,
				dna_nm: dnaConc,
			});
		} catch (e) {
			error = e.message;
		} finally {
			loading = false;
		}
	}

	function downloadCSV() {
		if (!results) return;
		const csv = exportCSV(results);
		const blob = new Blob([csv], { type: 'text/csv' });
		const url = URL.createObjectURL(blob);
		const a = document.createElement('a');
		a.href = url;
		a.download = 'primerdesignr_results.csv';
		a.click();
		URL.revokeObjectURL(url);
	}

	function copyIDTOrder() {
		if (!results) return;
		const lines = Object.entries(results.primers)
			.map(([name, p]) => `${name}\t${p.sequence}\t25nm\tSTD`);
		navigator.clipboard.writeText(lines.join('\n'));
	}

	/**
	 * Color-code a DNA sequence with spans per base
	 * @param {string} seq
	 * @returns {string} HTML string
	 */
	function colorSeq(seq) {
		return seq.split('').map(b => `<span class="${b}">${b}</span>`).join('');
	}

	/**
	 * Get the cross-dimer ΔG between two primers (order-insensitive)
	 */
	function getCrossDimer(a, b) {
		if (!results?.cross_dimers) return null;
		return results.cross_dimers.find(
			cd => (cd.primer_a === a && cd.primer_b === b) ||
			      (cd.primer_a === b && cd.primer_b === a)
		);
	}

	function dgClass(dg, threshold = -9) {
		if (dg < threshold) return 'status-danger';
		if (dg < threshold + 3) return 'status-warn';
		return 'status-ok';
	}
</script>

<svelte:head>
	<title>primerdesignr — Primer Analysis</title>
</svelte:head>

<!-- ── Input Section ──────────────────────────────────── -->
<section class="mb-6">
	<div class="flex items-end justify-between mb-3">
		<div>
			<h1 class="text-lg font-semibold">Primer Analysis</h1>
			<p class="text-xs text-[var(--color-text-muted)] mt-0.5">
				Paste sequences below. Accepts: NAME TAB SEQ, FASTA, or bare sequences.
			</p>
		</div>
		<span class="text-xs font-mono text-[var(--color-text-dim)]">
			{primerCount} primer{primerCount !== 1 ? 's' : ''} detected
		</span>
	</div>

	<textarea
		class="seq-input"
		rows="6"
		bind:value={input}
		placeholder={"F_primer\tATGCGATCGATCGATCGATCG\nR_primer\tCGATCGATCGATCGATCGAT\n\n...or paste FASTA, or bare sequences"}
	></textarea>

	<!-- Conditions -->
	<div class="flex items-center gap-6 mt-3 text-sm">
		<label class="flex items-center gap-2 text-[var(--color-text-muted)]">
			Na⁺
			<input type="number" bind:value={naConc} min="0" max="1000" step="10"
				class="w-16 bg-[var(--color-surface)] border border-[var(--color-surface-overlay)]
				       rounded px-2 py-1 text-[var(--color-text)] font-mono text-xs text-center
				       focus:border-[var(--color-info)] focus:outline-none" />
			<span class="text-xs text-[var(--color-text-dim)]">mM</span>
		</label>

		<label class="flex items-center gap-2 text-[var(--color-text-muted)]">
			Mg²⁺
			<input type="number" bind:value={mgConc} min="0" max="100" step="0.5"
				class="w-16 bg-[var(--color-surface)] border border-[var(--color-surface-overlay)]
				       rounded px-2 py-1 text-[var(--color-text)] font-mono text-xs text-center
				       focus:border-[var(--color-info)] focus:outline-none" />
			<span class="text-xs text-[var(--color-text-dim)]">mM</span>
		</label>

		<label class="flex items-center gap-2 text-[var(--color-text-muted)]">
			Primer
			<input type="number" bind:value={dnaConc} min="1" max="10000" step="50"
				class="w-20 bg-[var(--color-surface)] border border-[var(--color-surface-overlay)]
				       rounded px-2 py-1 text-[var(--color-text)] font-mono text-xs text-center
				       focus:border-[var(--color-info)] focus:outline-none" />
			<span class="text-xs text-[var(--color-text-dim)]">nM</span>
		</label>

		<div class="flex-1"></div>

		<button
			onclick={analyze}
			disabled={loading || primerCount === 0}
			class="px-4 py-1.5 rounded text-sm font-medium transition-all
			       bg-[var(--color-info)] text-[var(--color-bg)]
			       hover:brightness-110 disabled:opacity-40 disabled:cursor-not-allowed"
		>
			{loading ? 'Analyzing…' : 'Analyze'}
		</button>
	</div>

	{#if error}
		<div class="mt-3 px-3 py-2 rounded bg-red-500/10 border border-red-500/30 text-sm text-[var(--color-danger)]">
			{error}
		</div>
	{/if}
</section>

<!-- ── Results ────────────────────────────────────────── -->
{#if results}
	<!-- Summary bar -->
	<div class="mb-4 px-4 py-2.5 rounded-md bg-[var(--color-surface)] border border-[var(--color-surface-overlay)]
	            flex items-center justify-between">
		<div class="flex items-center gap-4 text-sm">
			<span class="{results.summary === 'All clear' ? 'status-ok' : 'status-warn'} font-medium">
				{results.summary}
			</span>
			<span class="text-[var(--color-text-dim)]">·</span>
			<span class="text-[var(--color-text-muted)]">
				Tm spread: <span class="font-mono">{results.tm_spread}°C</span>
			</span>
		</div>
		<div class="flex gap-2">
			<button onclick={downloadCSV}
				class="px-3 py-1 rounded text-xs font-medium
				       bg-[var(--color-surface-raised)] text-[var(--color-text-muted)]
				       hover:text-[var(--color-text)] transition-colors">
				Export CSV
			</button>
			<button onclick={copyIDTOrder}
				class="px-3 py-1 rounded text-xs font-medium
				       bg-[var(--color-surface-raised)] text-[var(--color-text-muted)]
				       hover:text-[var(--color-text)] transition-colors">
				Copy IDT Order
			</button>
		</div>
	</div>

	<!-- Primer table -->
	<div class="rounded-md border border-[var(--color-surface-overlay)] overflow-x-auto mb-6">
		<table class="data-table">
			<thead>
				<tr>
					<th></th>
					<th>Name</th>
					<th>Sequence</th>
					<th>Len</th>
					<th>GC</th>
					<th>Tm</th>
					<th>Hairpin (SL)</th>
					<th>Hairpin (MW)</th>
					<th>Self-dimer</th>
					<th>Structure</th>
				</tr>
			</thead>
			<tbody>
				{#each Object.entries(results.primers) as [name, p]}
					<tr>
						<td class="w-4">
							{#if p.warnings.length === 0}
								<span class="status-ok">✓</span>
							{:else}
								<span class="status-warn" title={p.warnings.join('\n')}>⚠</span>
							{/if}
						</td>
						<td class="font-medium text-sm whitespace-nowrap">{name}</td>
						<td class="seq whitespace-nowrap">{@html colorSeq(p.sequence)}</td>
						<td class="font-mono text-xs">{p.length}</td>
						<td class="font-mono text-xs">{(p.gc_content * 100).toFixed(0)}%</td>
						<td class="font-mono text-xs">{p.tm}°C</td>
						<td class="font-mono text-xs {dgClass(p.hairpin.dg_santalucia, -3)}">
							{p.hairpin.dg_santalucia}
						</td>
						<td class="font-mono text-xs {dgClass(p.hairpin.dg_mathews, -3)}">
							{p.hairpin.dg_mathews}
							{#if p.hairpin.engines_disagree}
								<span class="status-warn ml-1" title="SantaLucia/Mathews disagree — likely AT-closing stem">⚠</span>
							{/if}
						</td>
						<td class="font-mono text-xs {dgClass(p.homodimer.dg)}">
							{p.homodimer.dg}
						</td>
						<td class="dot-bracket whitespace-nowrap">{p.hairpin.dot_bracket}</td>
					</tr>

					<!-- Warning row -->
					{#if p.warnings.length > 0}
						<tr>
							<td></td>
							<td colspan="9" class="text-xs text-[var(--color-warn)] py-1">
								{p.warnings.join(' · ')}
							</td>
						</tr>
					{/if}
				{/each}
			</tbody>
		</table>
	</div>

	<!-- Cross-dimer matrix -->
	{#if results.cross_dimers.length > 0}
		<div class="mb-6">
			<h2 class="text-sm font-semibold mb-2 text-[var(--color-text-muted)] uppercase tracking-wider">
				Cross-dimer Matrix
				<span class="text-xs font-normal normal-case tracking-normal text-[var(--color-text-dim)]">
					ΔG kcal/mol · threshold: -9.0
				</span>
			</h2>
			<div class="rounded-md border border-[var(--color-surface-overlay)] overflow-x-auto">
				<table class="w-full border-collapse">
					<thead>
						<tr>
							<th class="matrix-cell text-[var(--color-text-dim)] text-xs"></th>
							{#each primerNames as name}
								<th class="matrix-cell text-[var(--color-text-muted)] text-xs font-medium">
									{name}
								</th>
							{/each}
						</tr>
					</thead>
					<tbody>
						{#each primerNames as row, i}
							<tr>
								<td class="matrix-cell text-[var(--color-text-muted)] text-xs font-medium text-left pl-3">
									{row}
								</td>
								{#each primerNames as col, j}
									{#if i === j}
										<td class="matrix-cell text-[var(--color-text-dim)]">—</td>
									{:else}
										{@const cd = getCrossDimer(row, col)}
										<td class="matrix-cell {cd ? dgClass(cd.dg) : ''}">
											{cd ? cd.dg.toFixed(1) : ''}
										</td>
									{/if}
								{/each}
							</tr>
						{/each}
					</tbody>
				</table>
			</div>
		</div>
	{/if}
{/if}
