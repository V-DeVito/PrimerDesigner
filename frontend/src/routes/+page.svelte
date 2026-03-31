<script>
	import { parsePrimerInput, analyzePrimers } from '$lib/api.js';

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
	let showConditions = $state(false);
	let expandedPrimer = $state(null);

	let parsedPrimers = $derived(parsePrimerInput(input));
	let primerCount = $derived(Object.keys(parsedPrimers).length);
	let primerNames = $derived(results ? Object.keys(results.primers) : []);

	let warningCount = $derived.by(() => {
		if (!results) return 0;
		return Object.values(results.primers).reduce((n, p) => n + (p.warnings.length > 0 ? 1 : 0), 0);
	});

	let tmRange = $derived.by(() => {
		if (!results) return null;
		const tms = Object.values(results.primers).map(p => p.tm);
		return { min: Math.min(...tms), max: Math.max(...tms) };
	});

	let gcRange = $derived.by(() => {
		if (!results) return null;
		const gcs = Object.values(results.primers).map(p => p.gc_content * 100);
		return { min: Math.min(...gcs), max: Math.max(...gcs) };
	});

	async function analyze() {
		const primers = parsedPrimers;
		if (Object.keys(primers).length === 0) { error = 'No valid primer sequences detected.'; return; }
		loading = true; error = ''; results = null; expandedPrimer = null;
		try { results = await analyzePrimers(primers, { na_mm: naConc, mg_mm: mgConc, dna_nm: dnaConc }); }
		catch (e) { error = e.message; }
		finally { loading = false; }
	}

	function colorSeq(seq) { return seq.split('').map(b => `<span class="${b}">${b}</span>`).join(''); }

	function getCrossDimer(a, b) {
		if (!results?.cross_dimers) return null;
		return results.cross_dimers.find(cd => (cd.primer_a === a && cd.primer_b === b) || (cd.primer_a === b && cd.primer_b === a));
	}

	// ── Color grading ──
	// Continuous severity: 0 (benign) → 1 (severe)
	function dgSeverity(dg, floor = -12) {
		if (dg === null || dg === undefined) return 0;
		const clamped = Math.max(floor, Math.min(0, dg));
		return Math.abs(clamped / floor);
	}

	// Background color: green → yellow → orange → red
	function dgBg(dg, floor = -12) {
		if (dg === null || dg === undefined) return '';
		const s = dgSeverity(dg, floor);
		let h, sat, light;
		if (s < 0.33) {
			h = 130 - (s / 0.33) * 80;
			sat = 55 + s * 30;
			light = 90 - s * 15;
		} else if (s < 0.66) {
			const t = (s - 0.33) / 0.33;
			h = 50 - t * 35;
			sat = 65 + t * 20;
			light = 85 - t * 10;
		} else {
			const t = (s - 0.66) / 0.34;
			h = 15 - t * 15;
			sat = 75 + t * 15;
			light = 75 - t * 10;
		}
		return `background: hsl(${h.toFixed(0)}, ${sat.toFixed(0)}%, ${light.toFixed(0)}%)`;
	}

	// Tm grading: distance from ideal 60°C
	function tmBg(tm) {
		const dist = Math.abs(tm - 60);
		if (dist > 8) return 'background: hsl(0, 75%, 85%)';
		if (dist > 4) return 'background: hsl(40, 70%, 88%)';
		return 'background: hsl(130, 50%, 90%)';
	}

	// GC% grading: distance from ideal 50%
	function gcBg(gc) {
		const dist = Math.abs(gc - 50);
		if (dist > 15) return 'background: hsl(0, 75%, 85%)';
		if (dist > 8) return 'background: hsl(40, 70%, 88%)';
		return 'background: hsl(130, 50%, 90%)';
	}

	function handleKeydown(e) {
		if ((e.metaKey || e.ctrlKey) && e.key === 'Enter') { e.preventDefault(); analyze(); }
	}
</script>

<svelte:head><title>Analyze — PrimerDesigner</title></svelte:head>
<svelte:window on:keydown={handleKeydown} />

<div class="pt-4 pb-12">

	<!-- ═══ INSTRUMENT PANEL ═══ -->
	<div class="border-[2.5px] border-[var(--color-ink)] enter"
		 style="box-shadow: var(--shadow-hard)">

		<!-- TITLE BAR -->
		<div class="flex items-center justify-between px-4 py-2 bg-[var(--color-ink)] text-[var(--color-cream)]">
			<span class="text-[13px] font-bold tracking-[0.08em] t-label">PRIMER ANALYSIS</span>
		</div>

		<!-- 3-PANEL GRID: input | summary | matrix -->
		<div class="grid grid-cols-1 lg:grid-cols-3">

			<!-- PANEL 1: INPUT -->
			<div class="border-b-[2.5px] lg:border-b-0 lg:border-r-[2.5px] border-[var(--color-ink)] flex flex-col">
				<textarea
					class="flex-1 w-full font-mono text-[13px] font-bold tracking-[0.04em] leading-[1.9]
						   bg-[var(--color-surface)] p-4 resize-none border-0 focus:outline-none min-h-[180px]"
					rows="7"
					bind:value={input}
					placeholder={"F_primer\tATGCGATCG…\nR_primer\tCGATCGATC…"}
					spellcheck="false"
					aria-label="Primer sequences"
				></textarea>
				<div class="flex items-center gap-3 px-4 py-2
							border-t-[2px] border-[var(--color-ink)] bg-[var(--color-surface-raised)]">
					<button onclick={analyze} disabled={loading || primerCount === 0}
							class="btn btn-primary !py-1 !px-4 !min-h-[32px] !text-[12px] !border-[2px] !shadow-none">
						{#if loading}
							<span class="spinner" style="width:12px;height:12px;border-width:2px"></span>
						{:else}
							Run
						{/if}
					</button>
					<button class="text-[12px] font-semibold text-[var(--color-text-muted)]
								   hover:text-[var(--color-text)] transition-colors"
							onclick={() => showConditions = !showConditions}>
						{showConditions ? '−' : '+'} Conditions
					</button>
				</div>
				{#if showConditions}
					<div class="flex items-center gap-3 px-4 py-2
								border-t border-[var(--color-border-subtle)] bg-[var(--color-surface-raised)] enter">
						<label class="flex items-center gap-1 text-[11px] text-[var(--color-text-secondary)]">
							<span class="font-semibold">Na+</span>
							<input type="number" bind:value={naConc} min="0" max="1000" step="10"
								class="w-12 text-[11px] text-center py-0.5 border border-[var(--color-border-subtle)]
									   bg-[var(--color-surface)] font-mono focus:outline-none focus:border-[var(--color-blue)]" />
						</label>
						<label class="flex items-center gap-1 text-[11px] text-[var(--color-text-secondary)]">
							<span class="font-semibold">Mg2+</span>
							<input type="number" bind:value={mgConc} min="0" max="100" step="0.5"
								class="w-12 text-[11px] text-center py-0.5 border border-[var(--color-border-subtle)]
									   bg-[var(--color-surface)] font-mono focus:outline-none focus:border-[var(--color-blue)]" />
						</label>
						<label class="flex items-center gap-1 text-[11px] text-[var(--color-text-secondary)]">
							<span class="font-semibold">DNA</span>
							<input type="number" bind:value={dnaConc} min="1" max="10000" step="50"
								class="w-14 text-[11px] text-center py-0.5 border border-[var(--color-border-subtle)]
									   bg-[var(--color-surface)] font-mono focus:outline-none focus:border-[var(--color-blue)]" />
							<span class="text-[var(--color-text-muted)]">nM</span>
						</label>
					</div>
				{/if}
				{#if error}
					<div class="px-4 py-2 bg-[var(--color-danger-bg)] text-[var(--color-danger)]
								text-[12px] font-semibold border-t border-[var(--color-border-subtle)]" role="alert">
						{error}
					</div>
				{/if}
			</div>

			<!-- PANEL 2: SUMMARY -->
			{#if results}
				<div class="border-b-[2.5px] lg:border-b-0 lg:border-r-[2.5px] border-[var(--color-ink)] flex flex-col">
					<!-- Verdict banner -->
					<div class="flex-1 flex items-center justify-center p-6
								{warningCount === 0 ? 'bg-[var(--color-ok)]' : 'bg-[var(--color-warn)]'}">
						<div class="text-white text-center">
							{#if warningCount === 0}
								<div class="text-4xl font-bold t-display leading-none">&#10003;</div>
								<div class="text-[11px] font-bold uppercase tracking-[0.15em] mt-2 opacity-90">All clear</div>
							{:else}
								<div class="text-4xl font-bold t-display leading-none">{warningCount}/{primerNames.length}</div>
								<div class="text-[11px] font-bold uppercase tracking-[0.15em] mt-2 opacity-90">
									primer{warningCount !== 1 ? 's' : ''} flagged
								</div>
							{/if}
						</div>
					</div>
					<!-- Metrics -->
					<div class="grid grid-cols-2 border-t-[2.5px] border-[var(--color-ink)]">
						<div class="px-4 py-3 border-r-[2.5px] border-[var(--color-ink)]"
							 style="{results.tm_spread > 8 ? 'background: hsl(0, 75%, 85%)' : results.tm_spread > 4 ? 'background: hsl(40, 70%, 88%)' : 'background: hsl(130, 50%, 90%)'}">
							<div class="text-[9px] font-bold uppercase tracking-[0.12em] text-[var(--color-text-muted)] mb-0.5">Tm spread</div>
							<div class="t-data text-xl font-bold leading-tight">{results.tm_spread.toFixed(1)}&deg;</div>
							{#if tmRange}
								<div class="t-data text-[11px] text-[var(--color-text-muted)] leading-tight">{tmRange.min.toFixed(1)}&ndash;{tmRange.max.toFixed(1)}&deg;C</div>
							{/if}
						</div>
						<div class="px-4 py-3"
							 style="{gcRange && (gcRange.min < 35 || gcRange.max > 65) ? 'background: hsl(0, 75%, 85%)' : gcRange && (gcRange.min < 40 || gcRange.max > 60) ? 'background: hsl(40, 70%, 88%)' : 'background: hsl(130, 50%, 90%)'}">
							<div class="text-[9px] font-bold uppercase tracking-[0.12em] text-[var(--color-text-muted)] mb-0.5">GC range</div>
							{#if gcRange}
								<div class="t-data text-xl font-bold leading-tight">{gcRange.min.toFixed(0)}&ndash;{gcRange.max.toFixed(0)}%</div>
								<div class="text-[11px] text-[var(--color-text-muted)] leading-tight">target 40&ndash;60%</div>
							{/if}
						</div>
					</div>
				</div>

				<!-- PANEL 3: DIMER MATRIX (self-dimer on diagonal) -->
				<div class="flex flex-col min-w-0">
					<div class="flex items-center justify-between px-3 py-1.5 bg-[var(--color-surface-raised)]
								border-b-[2px] border-[var(--color-ink)]">
						<span class="text-[10px] font-bold uppercase tracking-[0.1em] text-[var(--color-text-secondary)]">Dimer matrix &Delta;G</span>
						<span class="text-[9px] text-[var(--color-text-muted)]">kcal/mol</span>
					</div>
					<div class="flex-1 bg-[var(--color-surface)] p-2 overflow-x-auto flex items-center">
						<table class="w-full border-collapse" style="table-layout: fixed">
							<thead>
								<tr>
									<th class="w-[48px]"></th>
									{#each primerNames as name}
										<th class="text-[9px] font-bold uppercase tracking-wide text-[var(--color-text-muted)]
												   text-center whitespace-nowrap pb-1 px-0">{name}</th>
									{/each}
								</tr>
							</thead>
							<tbody>
								{#each primerNames as row, i}
									<tr>
										<td class="text-[9px] font-bold uppercase tracking-wide text-[var(--color-text-muted)]
												   pr-1.5 whitespace-nowrap text-right">{row}</td>
										{#each primerNames as col, j}
											{#if i === j}
												{@const selfDg = results.primers[row].homodimer.dg}
												<td class="border-[1.5px] border-[var(--color-ink)] text-[12px] text-center py-1 font-mono font-bold"
													style="{dgBg(selfDg)}"
													title="{row} self-dimer: {selfDg}">
													{selfDg}
												</td>
											{:else}
												{@const cd = getCrossDimer(row, col)}
												<td class="border-[1.5px] border-[var(--color-ink)] text-[12px] text-center py-1 font-mono font-bold"
													style="{cd ? dgBg(cd.dg) : ''}"
													title="{row} x {col}: {cd ? cd.dg.toFixed(2) : 'N/A'}">
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
			{:else if !loading}
				<div class="lg:col-span-2 flex items-center justify-center bg-[var(--color-surface-raised)] py-10">
					<div class="text-center">
						<div class="text-[13px] text-[var(--color-text-muted)] font-semibold">
							Paste primers and hit <span class="font-bold text-[var(--color-text-secondary)]">Run</span>
						</div>
						<div class="text-[11px] text-[var(--color-text-muted)] mt-1">
							Tm &middot; GC &middot; hairpin &middot; self-dimer &middot; cross-dimer
						</div>
					</div>
				</div>
			{/if}
		</div>
	</div>

	<!-- LOADING -->
	{#if loading}
		<div class="mt-4 enter">
			<div class="skeleton h-12 mb-1"></div>
			<div class="skeleton h-48"></div>
		</div>
	{/if}

	<!-- ═══ PRIMER TABLE ═══ -->
	{#if results}
		<div class="mt-4 border-[2.5px] border-[var(--color-ink)] enter" style="animation-delay: 60ms; box-shadow: var(--shadow-hard)">
			<!-- Table header -->
			<div class="grid px-4 py-2 bg-[var(--color-ink)] text-[var(--color-cream)]"
				 style="grid-template-columns: 1.2fr 3fr 1fr 1fr 1fr; gap: 0; align-items: center;">
				<span class="text-[10px] font-bold uppercase tracking-[0.1em]">Name</span>
				<span class="text-[10px] font-bold uppercase tracking-[0.1em]">Sequence</span>
				<span class="text-[10px] font-bold uppercase tracking-[0.1em] text-center">Tm</span>
				<span class="text-[10px] font-bold uppercase tracking-[0.1em] text-center">GC%</span>
				<span class="text-[10px] font-bold uppercase tracking-[0.1em] text-center">Hairpin</span>
			</div>

			{#each Object.entries(results.primers) as [name, p], idx}
				{@const hasWarnings = p.warnings.length > 0}
				{@const isExpanded = expandedPrimer === name}

				<div class="border-b border-[var(--color-border-subtle)]
							{hasWarnings ? 'bg-[var(--color-warn-bg)]' : 'bg-[var(--color-surface)]'}
							enter" style="animation-delay: {idx * 25}ms">
					<!-- Row -->
					<div class="grid px-4 py-2 cursor-pointer
								hover:bg-[var(--color-surface-raised)] transition-colors select-none"
						 style="grid-template-columns: 1.2fr 3fr 1fr 1fr 1fr; gap: 0; align-items: center;"
						 onclick={() => expandedPrimer = isExpanded ? null : name}
						 onkeydown={(e) => e.key === 'Enter' && (expandedPrimer = isExpanded ? null : name)}
						 role="button" tabindex="0" aria-expanded={isExpanded}>
						<span class="text-[13px] font-bold t-heading truncate">
							{#if hasWarnings}<span class="text-[var(--color-warn)] mr-1" title="{p.warnings.length} warning{p.warnings.length > 1 ? 's' : ''}">&#9888;</span>{/if}{name}
						</span>
						<span class="seq text-[12px] whitespace-nowrap overflow-hidden text-ellipsis min-w-0" title={p.sequence}>
							{@html colorSeq(p.sequence)}
						</span>
						<span class="t-data text-[13px] font-bold text-center px-1 py-0.5 mx-1"
							  style="{tmBg(p.tm)}">{p.tm}</span>
						<span class="t-data text-[13px] font-bold text-center px-1 py-0.5 mx-1"
							  style="{gcBg(p.gc_content * 100)}">{(p.gc_content * 100).toFixed(0)}</span>
						<span class="t-data text-[13px] font-bold text-center px-1 py-0.5 mx-1"
							  style="{dgBg(p.hairpin.dg_santalucia, -6)}">
							{p.hairpin.dg_santalucia}
						</span>
					</div>

					<!-- Expanded details: thermo blocks + warnings -->
					{#if isExpanded}
						<div class="px-4 pb-3 pt-1 enter border-t border-[var(--color-border-subtle)]">
							<div class="grid grid-cols-4 gap-0 mb-2 border border-[var(--color-border-subtle)]">
								<div class="px-3 py-2 border-r border-[var(--color-border-subtle)]"
									 style="{dgBg(p.hairpin.dg_santalucia, -6)}">
									<div class="text-[9px] font-bold uppercase tracking-[0.08em] text-[var(--color-text-muted)]">Hpin (SL)</div>
									<div class="t-data text-[14px] font-bold">{p.hairpin.dg_santalucia}</div>
								</div>
								<div class="px-3 py-2 border-r border-[var(--color-border-subtle)]"
									 style="{dgBg(p.hairpin.dg_mathews, -6)}">
									<div class="text-[9px] font-bold uppercase tracking-[0.08em] text-[var(--color-text-muted)]">Hpin (M)</div>
									<div class="t-data text-[14px] font-bold">
										{p.hairpin.dg_mathews}{#if p.hairpin.engines_disagree}<span class="text-[var(--color-warn)]"> &ne;</span>{/if}
									</div>
								</div>
								<div class="px-3 py-2 border-r border-[var(--color-border-subtle)]"
									 style="{dgBg(p.homodimer.dg)}">
									<div class="text-[9px] font-bold uppercase tracking-[0.08em] text-[var(--color-text-muted)]">Self-dim</div>
									<div class="t-data text-[14px] font-bold">{p.homodimer.dg}</div>
								</div>
								<div class="px-3 py-2"
									 style="{dgBg(p.homodimer.dg_3prime ?? 0)}">
									<div class="text-[9px] font-bold uppercase tracking-[0.08em] text-[var(--color-text-muted)]">3&prime; dim</div>
									<div class="t-data text-[14px] font-bold">
										{p.homodimer.dg_3prime ?? '—'}
									</div>
								</div>
							</div>
							{#if hasWarnings}
								{#each p.warnings as w}
									<p class="text-[11px] text-[var(--color-warn)] font-semibold leading-tight">&#9888; {w}</p>
								{/each}
							{/if}
						</div>
					{/if}
				</div>
			{/each}
		</div>
	{/if}
</div>
