<script>
	import { parsePrimerInput, analyzePrimers, exportCSV } from '$lib/api.js';

	// ── State ───────────────────────────────────────────
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
	let expandedPrimers = $state({});
	let toast = $state('');
	let toastKey = $state(0);

	// ── Derived ─────────────────────────────────────────
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

	// ── Actions ─────────────────────────────────────────
	async function analyze() {
		const primers = parsedPrimers;
		if (Object.keys(primers).length === 0) {
			error = 'No valid primer sequences detected.';
			return;
		}
		loading = true;
		error = '';
		results = null;
		expandedPrimers = {};
		try {
			results = await analyzePrimers(primers, { na_mm: naConc, mg_mm: mgConc, dna_nm: dnaConc });
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
		showToast('CSV downloaded');
	}

	function copyIDTOrder() {
		if (!results) return;
		const lines = Object.entries(results.primers)
			.map(([name, p]) => `${name}\t${p.sequence}\t25nm\tSTD`);
		navigator.clipboard.writeText(lines.join('\n'));
		showToast('Copied to clipboard');
	}

	function showToast(msg) {
		toast = msg;
		toastKey++;
		setTimeout(() => { toast = ''; }, 2200);
	}

	function togglePrimer(name) {
		expandedPrimers = { ...expandedPrimers, [name]: !expandedPrimers[name] };
	}

	function colorSeq(seq) {
		return seq.split('').map(b => `<span class="${b}">${b}</span>`).join('');
	}

	function getCrossDimer(a, b) {
		if (!results?.cross_dimers) return null;
		return results.cross_dimers.find(
			cd => (cd.primer_a === a && cd.primer_b === b) ||
			      (cd.primer_a === b && cd.primer_b === a)
		);
	}

	function dgColor(dg, threshold = -9) {
		if (dg < threshold) return 'val-danger';
		if (dg < threshold + 3) return 'val-warn';
		return 'val-ok';
	}

	function heatmapBg(dg) {
		if (dg === null || dg === undefined) return '';
		if (dg < -9) return 'background: rgba(239, 68, 68, 0.2)';
		if (dg < -6) return 'background: rgba(234, 179, 8, 0.15)';
		return 'background: rgba(34, 197, 94, 0.08)';
	}

	function handleKeydown(e) {
		if ((e.metaKey || e.ctrlKey) && e.key === 'Enter') {
			e.preventDefault();
			analyze();
		}
	}
</script>

<svelte:head>
	<title>Analyze — primerdesignr</title>
</svelte:head>

<svelte:window on:keydown={handleKeydown} />

<div class="max-w-7xl mx-auto px-5 py-6">
	<!-- ── Input Section ──────────────────────────────── -->
	<section class="mb-8">
		<!-- Header -->
		<div class="flex items-start justify-between mb-4">
			<div>
				<h1 class="text-xl font-semibold tracking-tight">Primer Analysis</h1>
				<p class="text-sm text-[var(--color-text-secondary)] mt-1">
					Paste sequences below — tab-separated, FASTA, or bare. &thinsp;
					<span class="text-[var(--color-text-muted)]">Up to 24 primers.</span>
				</p>
			</div>
			{#if primerCount > 0}
				<div class="flex items-center gap-1.5 bg-[var(--color-surface)] border border-[var(--color-border-subtle)]
				            rounded-full px-3 py-1">
					<div class="w-1.5 h-1.5 rounded-full bg-[var(--color-ok)]"></div>
					<span class="text-xs font-mono font-medium text-[var(--color-text-secondary)]">
						{primerCount} primer{primerCount !== 1 ? 's' : ''}
					</span>
				</div>
			{/if}
		</div>

		<!-- Textarea -->
		<textarea
			class="seq-input"
			rows={Math.min(Math.max(input.split('\n').length + 1, 4), 10)}
			bind:value={input}
			placeholder={"F_primer\tATGCGATCGATCGATCGATCG\nR_primer\tCGATCGATCGATCGATCGAT\n\n...or paste FASTA / bare sequences"}
			spellcheck="false"
		></textarea>

		<!-- Controls row -->
		<div class="flex items-center justify-between mt-3 gap-4">
			<div class="flex items-center gap-3">
				<!-- Conditions toggle -->
				<button
					class="btn btn-ghost btn-sm"
					onclick={() => showConditions = !showConditions}
				>
					<svg class="w-3.5 h-3.5 transition-transform {showConditions ? 'rotate-180' : ''}" fill="none" stroke="currentColor" viewBox="0 0 24 24">
						<path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M19 9l-7 7-7-7"/>
					</svg>
					Conditions
				</button>

				{#if showConditions}
					<div class="flex items-center gap-3 fade-in">
						<label class="flex items-center gap-1.5 text-xs text-[var(--color-text-secondary)]">
							Na&#8314;
							<input type="number" bind:value={naConc} min="0" max="1000" step="10"
								class="input input-mono w-14 text-xs text-center py-1 px-1.5" />
							<span class="text-[var(--color-text-muted)]">mM</span>
						</label>
						<label class="flex items-center gap-1.5 text-xs text-[var(--color-text-secondary)]">
							Mg&#178;&#8314;
							<input type="number" bind:value={mgConc} min="0" max="100" step="0.5"
								class="input input-mono w-14 text-xs text-center py-1 px-1.5" />
							<span class="text-[var(--color-text-muted)]">mM</span>
						</label>
						<label class="flex items-center gap-1.5 text-xs text-[var(--color-text-secondary)]">
							DNA
							<input type="number" bind:value={dnaConc} min="1" max="10000" step="50"
								class="input input-mono w-16 text-xs text-center py-1 px-1.5" />
							<span class="text-[var(--color-text-muted)]">nM</span>
						</label>
					</div>
				{/if}
			</div>

			<div class="flex items-center gap-2">
				<span class="text-[10px] text-[var(--color-text-muted)] hidden sm:inline">
					{#if !loading}Cmd+Enter{/if}
				</span>
				<button
					onclick={analyze}
					disabled={loading || primerCount === 0}
					class="btn btn-primary"
				>
					{#if loading}
						<span class="spinner"></span>
						Analyzing...
					{:else}
						<svg class="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
							<path stroke-linecap="round" stroke-linejoin="round" stroke-width="2"
							      d="M13 10V3L4 14h7v7l9-11h-7z"/>
						</svg>
						Analyze
					{/if}
				</button>
			</div>
		</div>

		{#if error}
			<div class="mt-3 px-4 py-2.5 rounded-lg bg-[var(--color-danger-soft)] border border-[var(--color-danger)]/20
			            text-sm text-[var(--color-danger)] flex items-center gap-2 fade-in">
				<svg class="w-4 h-4 shrink-0" fill="none" stroke="currentColor" viewBox="0 0 24 24">
					<path stroke-linecap="round" stroke-linejoin="round" stroke-width="2"
					      d="M12 8v4m0 4h.01M21 12a9 9 0 11-18 0 9 9 0 0118 0z"/>
				</svg>
				{error}
			</div>
		{/if}
	</section>

	<!-- ── Loading Skeleton ───────────────────────────── -->
	{#if loading}
		<section class="fade-in">
			<div class="grid gap-3">
				{#each Array(primerCount) as _}
					<div class="card">
						<div class="flex items-center gap-3 mb-3">
							<div class="skeleton w-20 h-4"></div>
							<div class="skeleton w-48 h-3"></div>
						</div>
						<div class="flex gap-6">
							<div class="skeleton w-12 h-8"></div>
							<div class="skeleton w-12 h-8"></div>
							<div class="skeleton w-12 h-8"></div>
							<div class="skeleton w-12 h-8"></div>
						</div>
					</div>
				{/each}
			</div>
		</section>
	{/if}

	<!-- ── Results ────────────────────────────────────── -->
	{#if results}
		<!-- Dashboard summary -->
		<section class="mb-6 fade-in">
			<div class="grid grid-cols-2 sm:grid-cols-4 gap-3">
				<!-- Status -->
				<div class="card flex items-center gap-3">
					{#if warningCount === 0}
						<div class="w-8 h-8 rounded-lg bg-[var(--color-ok-soft)] flex items-center justify-center">
							<svg class="w-4 h-4 text-[var(--color-ok)]" fill="none" stroke="currentColor" viewBox="0 0 24 24">
								<path stroke-linecap="round" stroke-linejoin="round" stroke-width="2.5" d="M5 13l4 4L19 7"/>
							</svg>
						</div>
						<div>
							<div class="text-sm font-semibold text-[var(--color-ok)]">All Clear</div>
							<div class="text-[11px] text-[var(--color-text-muted)]">No issues found</div>
						</div>
					{:else}
						<div class="w-8 h-8 rounded-lg bg-[var(--color-warn-soft)] flex items-center justify-center">
							<svg class="w-4 h-4 text-[var(--color-warn)]" fill="none" stroke="currentColor" viewBox="0 0 24 24">
								<path stroke-linecap="round" stroke-linejoin="round" stroke-width="2.5"
								      d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-2.5L13.732 4.5c-.77-.833-2.694-.833-3.464 0L3.34 16.5c-.77.833.192 2.5 1.732 2.5z"/>
							</svg>
						</div>
						<div>
							<div class="text-sm font-semibold text-[var(--color-warn)]">{warningCount} Warning{warningCount > 1 ? 's' : ''}</div>
							<div class="text-[11px] text-[var(--color-text-muted)]">Review flagged primers</div>
						</div>
					{/if}
				</div>

				<!-- Tm spread -->
				<div class="card">
					<div class="text-[11px] text-[var(--color-text-muted)] uppercase tracking-wider font-medium mb-1">Tm Spread</div>
					<div class="flex items-baseline gap-1.5">
						<span class="text-lg font-semibold font-mono {results.tm_spread > 5 ? 'text-[var(--color-warn)]' : 'text-[var(--color-text)]'}">
							{results.tm_spread.toFixed(1)}
						</span>
						<span class="text-xs text-[var(--color-text-muted)]">&#176;C</span>
					</div>
					{#if tmRange}
						<div class="text-[10px] text-[var(--color-text-muted)] font-mono mt-0.5">
							{tmRange.min.toFixed(1)} – {tmRange.max.toFixed(1)}&#176;C
						</div>
					{/if}
				</div>

				<!-- Primers count -->
				<div class="card">
					<div class="text-[11px] text-[var(--color-text-muted)] uppercase tracking-wider font-medium mb-1">Primers</div>
					<div class="text-lg font-semibold font-mono">{primerNames.length}</div>
					<div class="text-[10px] text-[var(--color-text-muted)] font-mono mt-0.5">
						{results.cross_dimers.length} cross-dimer pair{results.cross_dimers.length !== 1 ? 's' : ''}
					</div>
				</div>

				<!-- Export actions -->
				<div class="card flex flex-col justify-center gap-1.5">
					<button onclick={downloadCSV} class="btn btn-ghost btn-sm w-full">
						<svg class="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
							<path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M12 10v6m0 0l-3-3m3 3l3-3m2 8H7a2 2 0 01-2-2V5a2 2 0 012-2h5.586a1 1 0 01.707.293l5.414 5.414a1 1 0 01.293.707V19a2 2 0 01-2 2z"/>
						</svg>
						Export CSV
					</button>
					<button onclick={copyIDTOrder} class="btn btn-ghost btn-sm w-full">
						<svg class="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
							<path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M8 5H6a2 2 0 00-2 2v12a2 2 0 002 2h10a2 2 0 002-2v-1M8 5a2 2 0 002 2h2a2 2 0 002-2M8 5a2 2 0 012-2h2a2 2 0 012 2m0 0h2a2 2 0 012 2v3m2 4H10m0 0l3-3m-3 3l3 3"/>
						</svg>
						Copy IDT Order
					</button>
				</div>
			</div>
		</section>

		<!-- Primer Cards -->
		<section class="mb-8">
			<h2 class="text-xs font-semibold text-[var(--color-text-muted)] uppercase tracking-wider mb-3">
				Individual Primers
			</h2>
			<div class="grid gap-3">
				{#each Object.entries(results.primers) as [name, p], idx}
					{@const hasWarnings = p.warnings.length > 0}
					{@const isExpanded = expandedPrimers[name]}
					<div class="card {hasWarnings ? 'card-warn' : 'card-ok'} fade-in"
					     style="animation-delay: {idx * 50}ms">
						<!-- Card header -->
						<button class="w-full text-left" onclick={() => togglePrimer(name)}>
							<div class="flex items-center justify-between mb-2">
								<div class="flex items-center gap-2.5">
									<span class="text-sm font-semibold">{name}</span>
									{#if hasWarnings}
										<span class="badge badge-warn">{p.warnings.length} warning{p.warnings.length > 1 ? 's' : ''}</span>
									{:else}
										<span class="badge badge-ok">OK</span>
									{/if}
								</div>
								<svg class="w-4 h-4 text-[var(--color-text-muted)] transition-transform {isExpanded ? 'rotate-180' : ''}"
								     fill="none" stroke="currentColor" viewBox="0 0 24 24">
									<path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M19 9l-7 7-7-7"/>
								</svg>
							</div>

							<!-- Colored sequence -->
							<div class="seq overflow-x-auto pb-1 mb-3">{@html colorSeq(p.sequence)}</div>

							<!-- Metrics row -->
							<div class="flex flex-wrap gap-x-6 gap-y-2">
								<div class="metric">
									<span class="metric-label">Length</span>
									<span class="metric-value">{p.length} bp</span>
								</div>
								<div class="metric">
									<span class="metric-label">GC</span>
									<span class="metric-value">{(p.gc_content * 100).toFixed(0)}%</span>
								</div>
								<div class="metric">
									<span class="metric-label">Tm</span>
									<span class="metric-value">{p.tm}&#176;C</span>
								</div>
								<div class="metric">
									<span class="metric-label">Hairpin</span>
									<span class="metric-value {dgColor(p.hairpin.dg_santalucia, -3)}">
										{p.hairpin.dg_santalucia}
										{#if p.hairpin.engines_disagree}
											<span class="text-[var(--color-warn)] text-xs">*</span>
										{/if}
									</span>
								</div>
								<div class="metric">
									<span class="metric-label">Self-dimer</span>
									<span class="metric-value {dgColor(p.homodimer.dg)}">
										{p.homodimer.dg}
									</span>
								</div>
								{#if p.homodimer.has_3prime_risk}
									<div class="metric">
										<span class="metric-label">3' Risk</span>
										<span class="metric-value val-danger">
											{p.homodimer.dg_3prime}
										</span>
									</div>
								{/if}
							</div>
						</button>

						<!-- Expanded details -->
						{#if isExpanded}
							<div class="mt-4 pt-4 border-t border-[var(--color-border-subtle)] expand-enter">
								<!-- Detailed metrics grid -->
								<div class="grid grid-cols-2 sm:grid-cols-3 gap-4 mb-4">
									<div>
										<div class="text-[10px] uppercase tracking-wider text-[var(--color-text-muted)] mb-1 font-medium">
											Hairpin &#916;G (SantaLucia)
										</div>
										<div class="font-mono text-sm {dgColor(p.hairpin.dg_santalucia, -3)}">
											{p.hairpin.dg_santalucia} kcal/mol
										</div>
									</div>
									<div>
										<div class="text-[10px] uppercase tracking-wider text-[var(--color-text-muted)] mb-1 font-medium">
											Hairpin &#916;G (Mathews)
										</div>
										<div class="font-mono text-sm {dgColor(p.hairpin.dg_mathews, -3)}">
											{p.hairpin.dg_mathews} kcal/mol
											{#if p.hairpin.engines_disagree}
												<span class="badge badge-warn ml-1">engines disagree</span>
											{/if}
										</div>
									</div>
									<div>
										<div class="text-[10px] uppercase tracking-wider text-[var(--color-text-muted)] mb-1 font-medium">
											Self-dimer &#916;G
										</div>
										<div class="font-mono text-sm {dgColor(p.homodimer.dg)}">
											{p.homodimer.dg} kcal/mol
										</div>
									</div>
									{#if p.homodimer.dg_3prime !== undefined}
										<div>
											<div class="text-[10px] uppercase tracking-wider text-[var(--color-text-muted)] mb-1 font-medium">
												3' Self-dimer &#916;G
											</div>
											<div class="font-mono text-sm {p.homodimer.has_3prime_risk ? 'val-danger' : 'val-ok'}">
												{p.homodimer.dg_3prime} kcal/mol
											</div>
										</div>
									{/if}
								</div>

								<!-- Structure -->
								{#if p.hairpin.dot_bracket}
									<div class="mb-4">
										<div class="text-[10px] uppercase tracking-wider text-[var(--color-text-muted)] mb-1.5 font-medium">
											Secondary Structure
										</div>
										<div class="font-mono text-xs tracking-wider overflow-x-auto">
											<div class="seq text-[11px]">{@html colorSeq(p.sequence)}</div>
											<div class="dot-bracket">{p.hairpin.dot_bracket}</div>
										</div>
									</div>
								{/if}

								<!-- Warnings -->
								{#if hasWarnings}
									<div>
										<div class="text-[10px] uppercase tracking-wider text-[var(--color-text-muted)] mb-1.5 font-medium">
											Warnings
										</div>
										<div class="flex flex-col gap-1.5">
											{#each p.warnings as w}
												<div class="flex items-start gap-2 text-xs text-[var(--color-warn)]">
													<svg class="w-3.5 h-3.5 shrink-0 mt-0.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
														<path stroke-linecap="round" stroke-linejoin="round" stroke-width="2"
														      d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-2.5L13.732 4.5c-.77-.833-2.694-.833-3.464 0L3.34 16.5c-.77.833.192 2.5 1.732 2.5z"/>
													</svg>
													{w}
												</div>
											{/each}
										</div>
									</div>
								{/if}
							</div>
						{/if}
					</div>
				{/each}
			</div>
		</section>

		<!-- Cross-dimer Heatmap -->
		{#if results.cross_dimers.length > 0}
			<section class="mb-8 fade-in" style="animation-delay: 200ms">
				<div class="flex items-baseline justify-between mb-3">
					<h2 class="text-xs font-semibold text-[var(--color-text-muted)] uppercase tracking-wider">
						Cross-dimer Matrix
					</h2>
					<span class="text-[10px] text-[var(--color-text-muted)] font-mono">
						&#916;G kcal/mol &middot; threshold &#8804; -9.0
					</span>
				</div>
				<div class="card overflow-x-auto !p-4">
					<div class="inline-grid gap-1" style="grid-template-columns: max-content repeat({primerNames.length}, minmax(3.5rem, 1fr));">
						<!-- Header row -->
						<div></div>
						{#each primerNames as name}
							<div class="text-[10px] font-medium text-[var(--color-text-muted)] text-center px-1 whitespace-nowrap"
							     title={name}>
								{name}
							</div>
						{/each}

						<!-- Data rows -->
						{#each primerNames as row, i}
							<div class="text-[10px] font-medium text-[var(--color-text-muted)] flex items-center pr-3 whitespace-nowrap"
							     title={row}>
								{row}
							</div>
							{#each primerNames as col, j}
								{#if i === j}
									<div class="heatmap-cell text-[var(--color-text-muted)]/30 rounded-md">&#8212;</div>
								{:else if j > i}
									{@const cd = getCrossDimer(row, col)}
									<div class="heatmap-cell {cd ? dgColor(cd.dg) : ''} rounded-md cursor-default"
									     style="{cd ? heatmapBg(cd.dg) : ''}"
									     title="{row} x {col}: {cd ? cd.dg.toFixed(2) : 'N/A'} kcal/mol">
										{cd ? cd.dg.toFixed(1) : ''}
									</div>
								{:else}
									{@const cd = getCrossDimer(row, col)}
									<div class="heatmap-cell text-[var(--color-text-muted)]/40 rounded-md cursor-default"
									     title="{row} x {col}: {cd ? cd.dg.toFixed(2) : 'N/A'} kcal/mol">
										{cd ? cd.dg.toFixed(1) : ''}
									</div>
								{/if}
							{/each}
						{/each}
					</div>
				</div>
			</section>
		{/if}
	{/if}
</div>

<!-- Toast -->
{#if toast}
	{#key toastKey}
		<div class="toast">{toast}</div>
	{/key}
{/if}
