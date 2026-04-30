<script>
	import {
		analyzePrimers,
		checkGoldenGate,
		designPrimers,
		exportCSV,
		exportDesignCSV,
		normalizeSequence,
		parsePrimerInput
	} from '$lib/api.js';

	const sampleTemplate = `>example_template
TGACCTAAACGATGGTAGTAGGTTGGGAGCTTTCGAGAGGTCCGCCTTGGAAACGCGTTACTCGTGCGTGAATGTAGTGC
AAGAGGAGGGTCAGTCGTGCTGAGACTGGGACTCTAAAATTACTGAAGCTCTTCCCATCCTCTCTAAGGTTTCGTGAGAC
ACCACTTGGCGCTGGCCAGTACTGATCCCCTATTAGACCTATTGCCAAAACGTAGAAGGAACCTGCTCCAAAGCCCATGC
CATGTTTGCGGATAGAACCATCGGTTAAGCGCCCTGGGTCATTGGTGAACTGGGTGAAAGATTCAATCTGCTCTGGTGCC
CGGCCAGTCGCCATACGGAGCTCATTAGTATCGATTATAGAAGAGTTCTAGATCTTGCACTTGCGCGTCACAGAGTATAA
TTACTCCGTTACGGCATGGCGATGAAGCGATACTATAGTGAAATGAAACTTGTGGCCATCCAGGTCGTTACCGGCCGTGG
AACGGTGATATGTCAGACTTATAAGGTATATTGAGGGTTTAATAGCATTCCCGGCCAATTGTGTGTGCGATTGTTAAATG
GCAGCTTCCGGACTCACTAGGAACTATTGAACCTTATTGTCGGTGGGGATGTCCAATTTTAGAATTGAGCGTCGATGCAA
GGATCACAGTTTTAAAAGCAAGTTAAATCTAGGGTAAATAGCGGTGCTCTCATGCGTGTGGGTCGGGTGATGTTCAGAAA
ATTTGCCTCCGAGAATACACTCATGGGTAGGACTCGCACTACCTAAGATTCCGCGCGCAGCACCGTTCGAGATTCTGCCC`;

	let active = $state('design');

	let naConc = $state(50);
	let mgConc = $state(0);
	let dntpConc = $state(0);
	let dnaConc = $state(250);
	let showConditions = $state(false);

	let templateInput = $state(sampleTemplate);
	let designGoal = $state('exact');
	let productMin = $state(120);
	let productMax = $state(320);
	let primerCount = $state(5);
	let targetStart = $state('');
	let targetLength = $state('');
	let designLoading = $state(false);
	let designError = $state('');
	let designResults = $state(null);
	let expandedCandidate = $state(1);

	let analyzeInput = $state(`F_Chr\tGTCTTCACATCGGTTTGAAAGGAGG
R_Chr\tAACCCGCTCCGATTAAAGCTACTTT
F_TMM\tCCACAATGGAAGCAGATGTGATCGA
R_TMM\tACCACCCCAGTTATCCTGCTTTTCA
F_Kan\tTTTATTGATCTTGGGAGAAGCGGCA
R_Kan\tGTCCTGGGTTTCAAGCATTAGTCCA`);
	let analyzeLoading = $state(false);
	let analyzeError = $state('');
	let analyzeResults = $state(null);
	let expandedPrimer = $state(null);

	let overhangInput = $state('AACG AATG ATAG GCAA GCTG');
	let enzyme = $state('BsaI');
	let ggSequences = $state('');
	let ggLoading = $state(false);
	let ggError = $state('');
	let ggResults = $state(null);

	let templateLength = $derived(normalizeSequence(templateInput).length);
	let parsedPrimers = $derived(parsePrimerInput(analyzeInput));
	let analyzePrimerCount = $derived(Object.keys(parsedPrimers).length);
	let primerNames = $derived(analyzeResults ? Object.keys(analyzeResults.primers) : []);

	let warningCount = $derived.by(() => {
		if (!analyzeResults) return 0;
		return Object.values(analyzeResults.primers).reduce((n, p) => n + (p.warnings.length > 0 ? 1 : 0), 0);
	});

	let tmRange = $derived.by(() => {
		if (!analyzeResults) return null;
		const tms = Object.values(analyzeResults.primers).map((p) => p.tm);
		return { min: Math.min(...tms), max: Math.max(...tms) };
	});

	let tmSpread = $derived(analyzeResults ? analyzeResults.tm_spread : 0);

	function conditions() {
		return {
			na_mm: Number(naConc),
			mg_mm: Number(mgConc),
			dntp_mm: Number(dntpConc),
			dna_nm: Number(dnaConc)
		};
	}

	function parseOverhangs(raw) {
		return raw
			.split(/[\s,]+/)
			.map((s) => s.trim().toUpperCase())
			.filter(Boolean);
	}

	function parseSequenceMap(raw) {
		const parsed = parsePrimerInput(raw);
		return Object.keys(parsed).length ? parsed : null;
	}

	function optionalNumber(value) {
		return value === '' || value === null || value === undefined ? null : Number(value);
	}

	async function runDesign() {
		designLoading = true;
		designError = '';
		designResults = null;
		expandedCandidate = 1;

		if (designGoal === 'targeted' && (targetStart === '' || targetLength === '')) {
			designError = 'Required region coordinates are needed.';
			designLoading = false;
			return;
		}

		try {
			designResults = await designPrimers({
				template: templateInput,
				design_mode: designGoal,
				product_min: Number(productMin),
				product_max: Number(productMax),
				primer_count: designGoal === 'exact' ? 1 : Number(primerCount),
				target_start: designGoal === 'exact' ? null : optionalNumber(targetStart),
				target_length: designGoal === 'exact' ? null : optionalNumber(targetLength),
				conditions: conditions()
			});
		} catch (e) {
			designError = e.message;
		} finally {
			designLoading = false;
		}
	}

	async function runAnalyze() {
		const primers = parsedPrimers;
		if (Object.keys(primers).length === 0) {
			analyzeError = 'No valid primer sequences detected.';
			return;
		}

		analyzeLoading = true;
		analyzeError = '';
		analyzeResults = null;
		expandedPrimer = null;

		try {
			analyzeResults = await analyzePrimers(primers, conditions());
		} catch (e) {
			analyzeError = e.message;
		} finally {
			analyzeLoading = false;
		}
	}

	async function runGoldenGate() {
		ggLoading = true;
		ggError = '';
		ggResults = null;

		try {
			ggResults = await checkGoldenGate(parseOverhangs(overhangInput), enzyme, parseSequenceMap(ggSequences));
		} catch (e) {
			ggError = e.message;
		} finally {
			ggLoading = false;
		}
	}

	function getCrossDimer(a, b) {
		if (!analyzeResults?.cross_dimers) return null;
		return analyzeResults.cross_dimers.find(
			(cd) => (cd.primer_a === a && cd.primer_b === b) || (cd.primer_a === b && cd.primer_b === a)
		);
	}

	// Cross-dimer heatmap: more negative ΔG → stronger blush wash.
	function heatBg(value, floor = -6, ceiling = -1) {
		if (value === null || value === undefined) return '';
		const t = Math.max(0, Math.min(1, (value - floor) / (ceiling - floor)));
		const opacity = ((1 - t) * 0.85 * 0.55).toFixed(3);
		return `background: rgba(198, 110, 110, ${opacity})`;
	}

	function downloadCsv(filename, csv) {
		const blob = new Blob([csv], { type: 'text/csv;charset=utf-8' });
		const url = URL.createObjectURL(blob);
		const a = document.createElement('a');
		a.href = url;
		a.download = filename;
		a.click();
		URL.revokeObjectURL(url);
	}

	async function copyPrimerPair(candidate) {
		const text = [
			`Forward\t${candidate.forward.sequence}`,
			`Reverse\t${candidate.reverse.sequence}`
		].join('\n');
		await navigator.clipboard?.writeText(text);
	}
</script>

<svelte:head><title>Primer Workbench</title></svelte:head>

<!-- Tabs -->
<div class="flex" style="gap: 28px; border-bottom: 1px solid var(--line); margin-bottom: var(--gap-section);">
	<button class="pw-tab" data-active={active === 'design'} onclick={() => active = 'design'}>Design</button>
	<button class="pw-tab" data-active={active === 'analyze'} onclick={() => active = 'analyze'}>Analyze</button>
	<button class="pw-tab" data-active={active === 'golden'} onclick={() => active = 'golden'}>Golden Gate</button>
</div>

<!-- ─────────── DESIGN ─────────── -->
{#if active === 'design'}
	<div class="design-grid items-start">
		<!-- Template editor -->
		<div class="pw-card">
			<div class="flex items-center justify-between" style="margin-bottom: 12px;">
				<span class="pw-section-label">Template</span>
				<span class="pw-num" style="font-size: 12px; color: var(--muted);">{templateLength} bp</span>
			</div>
			<textarea
				class="pw-textarea"
				style="min-height: 280px;"
				bind:value={templateInput}
				spellcheck="false"
			></textarea>
		</div>

		<!-- Constraints -->
		<div class="pw-card">
			<div style="margin-bottom: 14px;">
				<span class="pw-section-label">Constraints</span>
			</div>

			<div style="margin-bottom: 14px;">
				<div class="pw-field-label" style="margin-bottom: 6px;">Design goal</div>
				<div class="pw-select-wrap">
					<select class="pw-select" bind:value={designGoal} style="padding-right: 28px;">
						<option value="exact">Exact submitted bounds</option>
						<option value="targeted">Required region</option>
						<option value="amplicon">Exploratory internal amplicon</option>
					</select>
				</div>
			</div>

			<div class="grid" style="grid-template-columns: 1fr 1fr; gap: 12px; margin-bottom: 14px;">
				<div>
					<div class="pw-field-label" style="margin-bottom: 6px;">Min product</div>
					<input class="pw-input" type="number" min="40" max="5000" bind:value={productMin} disabled={designGoal === 'exact'} />
				</div>
				<div>
					<div class="pw-field-label" style="margin-bottom: 6px;">Max product</div>
					<input class="pw-input" type="number" min="40" max="5000" bind:value={productMax} disabled={designGoal === 'exact'} />
				</div>
			</div>

			<div style="margin-bottom: 14px;">
				<div class="pw-field-label" style="margin-bottom: 6px;">Candidates</div>
				<input class="pw-input" type="number" min="1" max="20"
					value={designGoal === 'exact' ? 1 : primerCount}
					disabled={designGoal === 'exact'}
					oninput={(e) => primerCount = e.currentTarget.value} />
			</div>

			{#if designGoal !== 'exact'}
				<div class="grid" style="grid-template-columns: 1fr 1fr; gap: 12px; margin-bottom: 18px;">
					<div>
						<div class="pw-field-label" style="margin-bottom: 6px;">{designGoal === 'targeted' ? 'Req. start' : 'Target start'}</div>
						<input class="pw-input" type="number" min="0" bind:value={targetStart} />
					</div>
					<div>
						<div class="pw-field-label" style="margin-bottom: 6px;">{designGoal === 'targeted' ? 'Req. length' : 'Target length'}</div>
						<input class="pw-input" type="number" min="1" bind:value={targetLength} />
					</div>
				</div>
			{/if}

			<button class="pw-btn-primary" style="width: 100%; margin-top: 4px;"
				disabled={designLoading || templateLength < 40}
				onclick={runDesign}>
				{designLoading ? 'Designing…' : 'Design primers'}
			</button>
			{#if designError}
				<p style="margin-top: 12px; font-size: 12px; color: var(--blush-fg);">{designError}</p>
			{/if}
		</div>
	</div>

	{#if designResults}
		<div style="margin-top: var(--gap-section);">
			<div class="pw-card pw-card-flush">
				<div class="flex items-center justify-between" style="padding: var(--pad-card) var(--pad-card) 12px;">
					<span class="pw-section-label">Candidates</span>
					<button class="pw-btn-ghost"
						onclick={() => downloadCsv('primerdesign-candidates.csv', exportDesignCSV(designResults))}>
						Export CSV
					</button>
				</div>

				{#if designResults.warnings.length || designResults.candidates.length === 0}
					<div style="margin: 0 var(--pad-card) 12px; padding: 12px 14px; background: var(--surface2); border: 1px solid var(--line-soft); border-radius: 10px;">
						{#if designResults.warnings.length}
							{#each designResults.warnings as warning}
								<p style="font-size: 12.5px; color: var(--ink-soft); margin: 2px 0;">{warning}</p>
							{/each}
						{/if}
						{#if designResults.candidates.length === 0 && Object.keys(designResults.primer3_explain).length}
							<div class="grid" style="margin-top: 10px; gap: 6px; font-size: 11.5px; color: var(--muted); grid-template-columns: 1fr 1fr;">
								{#each Object.entries(designResults.primer3_explain) as [key, value]}
									<p><span class="pw-num" style="color: var(--muted-soft);">{key}</span>: {value}</p>
								{/each}
							</div>
						{/if}
					</div>
				{/if}

				{#if designResults.candidates.length}
					<table style="width: 100%; border-collapse: collapse;">
						<thead>
							<tr>
								{#each ['Rank', 'Forward', 'Reverse', 'Product', 'Tm Δ', 'Dimer'] as h, i}
									<th style="padding: 14px 16px; text-align: {i >= 3 ? 'right' : 'left'}; border-bottom: 1px solid var(--line);"
										class="pw-eyebrow">{h}</th>
								{/each}
							</tr>
						</thead>
						<tbody>
							{#each designResults.candidates as c, i}
								<tr style="background: {i % 2 === 1 ? 'var(--row-alt)' : 'transparent'}; cursor: pointer; border-bottom: 1px solid var(--line-soft);"
									onclick={() => expandedCandidate = expandedCandidate === c.rank ? null : c.rank}>
									<td class="pw-num" style="padding: 13px 16px; font-size: 13.5px; color: var(--ink);">{c.rank}</td>
									<td class="pw-num" style="padding: 13px 16px; font-size: 13px; color: var(--ink-soft); white-space: nowrap; overflow: hidden; text-overflow: ellipsis; max-width: 0;">{c.forward.sequence}</td>
									<td class="pw-num" style="padding: 13px 16px; font-size: 13px; color: var(--ink-soft); white-space: nowrap; overflow: hidden; text-overflow: ellipsis; max-width: 0;">{c.reverse.sequence}</td>
									<td class="pw-num" style="padding: 13px 16px; font-size: 13.5px; text-align: right; color: var(--ink);">{c.product_size}</td>
									<td class="pw-num" style="padding: 13px 16px; font-size: 13.5px; text-align: right; color: var(--ink);">{c.tm_difference.toFixed(1)}</td>
									<td class="pw-num" style="padding: 13px 16px; font-size: 13.5px; text-align: right; color: var(--ink); {heatBg(c.heterodimer.dg)}">{c.heterodimer.dg.toFixed(2)}</td>
								</tr>
								{#if expandedCandidate === c.rank}
									<tr>
										<td colspan="6" style="padding: 16px 20px; background: var(--surface2); border-bottom: 1px solid var(--line-soft);">
											<div class="grid" style="grid-template-columns: 1fr 1fr 1fr; gap: 18px;">
												<div>
													<div class="pw-field-label" style="margin-bottom: 6px;">Forward</div>
													<p class="pw-num" style="font-size: 12.5px; color: var(--ink);">{c.forward_coords.start}-{c.forward_coords.end} ({c.forward_coords.strand})</p>
													<p style="font-size: 12px; color: var(--muted); margin-top: 4px;">Tm {c.forward.tm} °C · GC {(c.forward.gc_content * 100).toFixed(0)}%</p>
												</div>
												<div>
													<div class="pw-field-label" style="margin-bottom: 6px;">Reverse</div>
													<p class="pw-num" style="font-size: 12.5px; color: var(--ink);">{c.reverse_coords.start}-{c.reverse_coords.end} ({c.reverse_coords.strand})</p>
													<p style="font-size: 12px; color: var(--muted); margin-top: 4px;">Tm {c.reverse.tm} °C · GC {(c.reverse.gc_content * 100).toFixed(0)}%</p>
												</div>
												<div>
													<div class="pw-field-label" style="margin-bottom: 6px;">Actions</div>
													<button class="pw-btn-ghost" onclick={(e) => { e.stopPropagation(); copyPrimerPair(c); }}>Copy pair</button>
												</div>
											</div>
											{#if c.warnings.length}
												<div style="margin-top: 14px; padding-top: 12px; border-top: 1px solid var(--line-soft);">
													{#each c.warnings as warning}
														<p style="display: flex; align-items: center; gap: 8px; font-size: 12.5px; color: var(--ink-soft); margin: 2px 0;">
															<span class="pw-dot pw-dot-warn"></span>{warning}
														</p>
													{/each}
												</div>
											{/if}
										</td>
									</tr>
								{/if}
							{/each}
						</tbody>
					</table>
				{/if}
			</div>
		</div>
	{/if}
{/if}

<!-- ─────────── ANALYZE ─────────── -->
{#if active === 'analyze'}
	<div style="display: flex; flex-direction: column; gap: var(--gap-section);">
		<!-- Primer set editor -->
		<div class="pw-card">
			<div style="margin-bottom: 12px;">
				<span class="pw-section-label">Primer set</span>
			</div>
			<textarea
				class="pw-textarea"
				style="min-height: 160px; line-height: 1.85;"
				bind:value={analyzeInput}
				spellcheck="false"
				placeholder={"F_primer\tATGCGATCG…\nR_primer\tCGATCGATC…"}
			></textarea>

			<div class="flex items-center" style="gap: 10px; margin-top: 14px;">
				<button class="pw-btn-primary"
					disabled={analyzeLoading || analyzePrimerCount === 0}
					onclick={runAnalyze}>
					{analyzeLoading ? 'Analyzing…' : 'Analyze set'}
				</button>
				{#if analyzeResults}
					<button class="pw-btn-ghost"
						onclick={() => downloadCsv('primerdesign-analysis.csv', exportCSV(analyzeResults))}>
						Export CSV
					</button>
				{/if}
				<button class="pw-btn-quiet" style="margin-left: auto;"
					onclick={() => showConditions = !showConditions}>
					{showConditions ? 'Hide conditions' : 'Conditions'}
				</button>
			</div>

			{#if showConditions}
				<div class="flex flex-wrap" style="gap: 18px; padding: 12px 14px; border-radius: 10px; background: var(--surface2); border: 1px solid var(--line-soft); margin-top: 14px; font-family: var(--font-mono); font-size: 12.5px; color: var(--muted); letter-spacing: 0.02em;">
					<label class="flex items-center" style="gap: 6px;">
						<span style="color: var(--ink);">Na+</span>
						<input class="pw-input" type="number" min="0" max="1000" step="10" bind:value={naConc}
							style="width: 70px; height: 26px; padding: 0 6px; font-size: 12px;" />
						<span>mM</span>
					</label>
					<label class="flex items-center" style="gap: 6px;">
						<span style="color: var(--ink);">Mg2+</span>
						<input class="pw-input" type="number" min="0" max="100" step="0.5" bind:value={mgConc}
							style="width: 70px; height: 26px; padding: 0 6px; font-size: 12px;" />
						<span>mM</span>
					</label>
					<label class="flex items-center" style="gap: 6px;">
						<span style="color: var(--ink);">dNTP</span>
						<input class="pw-input" type="number" min="0" max="20" step="0.1" bind:value={dntpConc}
							style="width: 70px; height: 26px; padding: 0 6px; font-size: 12px;" />
						<span>mM</span>
					</label>
					<label class="flex items-center" style="gap: 6px;">
						<span style="color: var(--ink);">Primer</span>
						<input class="pw-input" type="number" min="1" max="10000" step="50" bind:value={dnaConc}
							style="width: 80px; height: 26px; padding: 0 6px; font-size: 12px;" />
						<span>nM</span>
					</label>
				</div>
			{/if}
			{#if analyzeError}
				<p style="margin-top: 12px; font-size: 12px; color: var(--blush-fg);">{analyzeError}</p>
			{/if}
		</div>

		{#if analyzeResults}
			<!-- Metric cards -->
			<div class="grid" style="grid-template-columns: repeat(3, 1fr); gap: var(--gap-col);">
				<div class="pw-card">
					<div class="flex items-center justify-between" style="margin-bottom: 12px;">
						<span class="pw-eyebrow">Warnings</span>
						<span class="pw-dot {warningCount > 0 ? 'pw-dot-warn' : 'pw-dot-ok'}"></span>
					</div>
					<div class="flex items-baseline" style="gap: 4px;">
						<span class="pw-num" style="font-size: 29px; color: var(--ink); letter-spacing: -0.01em;">{warningCount}</span>
						<span class="pw-num" style="font-size: 15px; color: var(--muted);">/{primerNames.length}</span>
					</div>
				</div>
				<div class="pw-card">
					<div class="flex items-center justify-between" style="margin-bottom: 12px;">
						<span class="pw-eyebrow">Tm spread</span>
						<span class="pw-dot {tmSpread > 3 ? 'pw-dot-warn' : 'pw-dot-ok'}"></span>
					</div>
					<div class="flex items-baseline" style="gap: 4px;">
						<span class="pw-num" style="font-size: 29px; color: var(--ink); letter-spacing: -0.01em;">{tmSpread.toFixed(1)}</span>
						<span class="pw-num" style="font-size: 15px; color: var(--muted);">°C</span>
					</div>
				</div>
				<div class="pw-card">
					<div style="margin-bottom: 12px;">
						<span class="pw-eyebrow">Range</span>
					</div>
					<div class="flex items-baseline" style="gap: 4px;">
						<span class="pw-num" style="font-size: 22px; color: var(--ink); letter-spacing: -0.01em;">
							{tmRange ? `${tmRange.min.toFixed(1)}–${tmRange.max.toFixed(1)}` : '–'}
						</span>
						<span class="pw-num" style="font-size: 15px; color: var(--muted);">°C</span>
					</div>
				</div>
			</div>

			<!-- Primer table -->
			<div class="pw-card pw-card-flush">
				<table style="width: 100%; border-collapse: collapse; table-layout: fixed;">
					<colgroup>
						<col style="width: 14%;" />
						<col style="width: 35%;" />
						<col style="width: 11%;" />
						<col style="width: 11%;" />
						<col style="width: 13%;" />
						<col style="width: 16%;" />
					</colgroup>
					<thead>
						<tr>
							{#each [['Name', 'left'], ['Sequence', 'left'], ['Tm', 'right'], ['GC', 'right'], ['Hairpin', 'right'], ['Self-dimer', 'right']] as [h, align]}
								<th class="pw-eyebrow" style="padding: 14px 16px; text-align: {align}; border-bottom: 1px solid var(--line);">{h}</th>
							{/each}
						</tr>
					</thead>
					<tbody>
						{#each Object.entries(analyzeResults.primers) as [name, p], i}
							{@const isExpanded = expandedPrimer === name}
							<tr style="background: {i % 2 === 1 ? 'var(--row-alt)' : 'transparent'}; cursor: pointer; border-bottom: 1px solid var(--line-soft);"
								onclick={() => expandedPrimer = isExpanded ? null : name}>
								<td style="padding: 13px 16px; font-size: 13.5px; color: var(--ink); white-space: nowrap;">
									<span class="inline-flex items-center" style="gap: 8px;">
										{#if p.warnings.length}<span class="pw-dot pw-dot-warn"></span>{/if}
										<span class="pw-num">{name}</span>
									</span>
								</td>
								<td class="pw-num" style="padding: 13px 16px; font-size: 13px; color: var(--ink-soft); white-space: nowrap; overflow: hidden; text-overflow: ellipsis;">{p.sequence}</td>
								<td class="pw-num" style="padding: 13px 16px; font-size: 13.5px; text-align: right; color: var(--ink);">{p.tm.toFixed(1)}</td>
								<td class="pw-num" style="padding: 13px 16px; font-size: 13.5px; text-align: right; color: var(--ink);">{(p.gc_content * 100).toFixed(0)}%</td>
								<td class="pw-num" style="padding: 13px 16px; font-size: 13.5px; text-align: right; color: var(--ink); {heatBg(p.hairpin.dg_santalucia, -6, -1)}">{p.hairpin.dg_santalucia}</td>
								<td class="pw-num" style="padding: 13px 16px; font-size: 13.5px; text-align: right; color: var(--ink); {heatBg(p.homodimer.dg)}">{p.homodimer.dg}</td>
							</tr>
							{#if isExpanded}
								<tr>
									<td colspan="6" style="padding: 14px 20px; background: var(--surface2); border-bottom: 1px solid var(--line-soft);">
										<div style="font-size: 12.5px;">
											<span style="color: var(--muted);">Mathews hairpin</span>
											<span class="pw-num" style="margin-left: 8px; color: var(--ink);">{p.hairpin.dg_mathews}</span>
										</div>
										{#if p.warnings.length}
											<div style="margin-top: 12px; padding-top: 10px; border-top: 1px solid var(--line-soft);">
												{#each p.warnings as warning}
													<p style="display: flex; align-items: center; gap: 8px; font-size: 12.5px; color: var(--ink-soft); margin: 2px 0;">
														<span class="pw-dot pw-dot-warn"></span>{warning}
													</p>
												{/each}
											</div>
										{/if}
									</td>
								</tr>
							{/if}
						{/each}
					</tbody>
				</table>
			</div>

			<!-- Cross-dimer matrix -->
			<div class="pw-card pw-card-flush">
				<div style="padding: var(--pad-card) var(--pad-card) 12px;">
					<span class="pw-section-label">Cross-dimer matrix</span>
				</div>
				<div style="overflow-x: auto;">
					<table style="width: 100%; border-collapse: collapse; table-layout: fixed;">
						<thead>
							<tr>
								<th class="pw-eyebrow" style="width: 96px; padding: 12px;"></th>
								{#each primerNames as name}
									<th class="pw-eyebrow" style="padding: 12px; text-align: center;">{name}</th>
								{/each}
							</tr>
						</thead>
						<tbody>
							{#each primerNames as row, i}
								<tr>
									<td class="pw-eyebrow" style="padding: 12px; text-align: right; border-right: 1px solid var(--line-soft); border-top: 1px solid var(--line-soft);">{row}</td>
									{#each primerNames as col, j}
										{#if i === j}
											{@const selfDg = analyzeResults.primers[row].homodimer.dg}
											<td class="pw-num"
												style="padding: 0; text-align: center; font-size: 13.5px; border-top: 1px solid var(--line-soft); {j < primerNames.length - 1 ? 'border-right: 1px solid var(--line-soft);' : ''} {heatBg(selfDg)}">
												<div style="padding: 13px 6px; font-style: italic; color: var(--muted);">{selfDg}</div>
											</td>
										{:else}
											{@const cd = getCrossDimer(row, col)}
											<td class="pw-num"
												style="padding: 0; text-align: center; font-size: 13.5px; border-top: 1px solid var(--line-soft); {j < primerNames.length - 1 ? 'border-right: 1px solid var(--line-soft);' : ''} {cd ? heatBg(cd.dg) : ''}">
												<div style="padding: 13px 6px; color: var(--ink);">{cd ? cd.dg.toFixed(2) : ''}</div>
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
	</div>
{/if}

<!-- ─────────── GOLDEN GATE ─────────── -->
{#if active === 'golden'}
	<div style="display: flex; flex-direction: column; gap: var(--gap-section);">
		<div class="pw-card">
			<div class="grid" style="grid-template-columns: 1fr 1fr 1fr; gap: var(--gap-col);">
				<div>
					<div class="pw-field-label" style="margin-bottom: 6px;">Type-IIS enzyme</div>
					<div class="pw-select-wrap">
						<select class="pw-select" bind:value={enzyme} style="padding-right: 28px;">
							<option>BsaI</option>
							<option>BsmBI</option>
							<option>BbsI</option>
							<option>SapI</option>
							<option>BpiI</option>
						</select>
					</div>
				</div>
				<div>
					<div class="pw-field-label" style="margin-bottom: 6px;">Overhang length</div>
					<div class="pw-select-wrap">
						<select class="pw-select" disabled style="padding-right: 28px;">
							<option>4 nt</option>
						</select>
					</div>
				</div>
				<div>
					<div class="pw-field-label" style="margin-bottom: 6px;">Strategy</div>
					<div class="pw-select-wrap">
						<select class="pw-select" disabled style="padding-right: 28px;">
							<option>One-pot</option>
						</select>
					</div>
				</div>
			</div>
		</div>

		<div class="pw-card">
			<div style="margin-bottom: 12px;">
				<span class="pw-section-label">Overhangs</span>
			</div>
			<textarea class="pw-textarea" style="min-height: 100px;" bind:value={overhangInput} spellcheck="false"></textarea>

			<div style="margin-top: 14px; margin-bottom: 12px;">
				<span class="pw-section-label">Sequences to scan</span>
			</div>
			<textarea class="pw-textarea" style="min-height: 100px;" bind:value={ggSequences}
				placeholder={">insert\nGGTCTC…"} spellcheck="false"></textarea>

			<div class="flex" style="gap: 10px; margin-top: 14px;">
				<button class="pw-btn-primary"
					disabled={ggLoading || parseOverhangs(overhangInput).length < 2}
					onclick={runGoldenGate}>
					{ggLoading ? 'Checking…' : 'Check overhangs'}
				</button>
			</div>
			{#if ggError}
				<p style="margin-top: 12px; font-size: 12px; color: var(--blush-fg);">{ggError}</p>
			{/if}
		</div>

		{#if ggResults}
			<div class="pw-card pw-card-flush">
				<div style="padding: var(--pad-card) var(--pad-card) 12px;">
					<span class="pw-section-label">Report</span>
				</div>
				<table style="width: 100%; border-collapse: collapse;">
					<thead>
						<tr>
							{#each [['Overhang', 'left'], ['Reverse complement', 'left'], ['GC', 'right'], ['Fidelity set', 'left'], ['Warnings', 'left']] as [h, align]}
								<th class="pw-eyebrow" style="padding: 14px 16px; text-align: {align}; border-bottom: 1px solid var(--line);">{h}</th>
							{/each}
						</tr>
					</thead>
					<tbody>
						{#each ggResults.overhangs as oh, i}
							<tr style="background: {i % 2 === 1 ? 'var(--row-alt)' : 'transparent'}; border-bottom: 1px solid var(--line-soft);">
								<td style="padding: 13px 16px;">
									<span class="pw-pill-sky">{oh.overhang}</span>
								</td>
								<td class="pw-num" style="padding: 13px 16px; font-size: 13.5px; color: var(--ink-soft);">{oh.reverse_complement}</td>
								<td class="pw-num" style="padding: 13px 16px; font-size: 13.5px; text-align: right; color: var(--ink);">{(oh.gc_content * 100).toFixed(0)}%</td>
								<td style="padding: 13px 16px; font-size: 13.5px; color: var(--ink-soft);">{oh.is_in_fidelity_set ? 'NEB high-fidelity' : 'custom'}</td>
								<td style="padding: 13px 16px; font-size: 12.5px; color: var(--ink-soft);">
									{#if oh.warnings.length}
										<span class="inline-flex items-center" style="gap: 6px;">
											<span class="pw-dot pw-dot-warn"></span>{oh.warnings.join('; ')}
										</span>
									{:else}
										<span style="color: var(--muted);">—</span>
									{/if}
								</td>
							</tr>
						{/each}
					</tbody>
				</table>
				{#if ggResults.warnings.length}
					<div style="margin: 12px var(--pad-card) var(--pad-card); padding: 12px 14px; background: var(--surface2); border: 1px solid var(--line-soft); border-radius: 10px;">
						<div class="pw-eyebrow" style="margin-bottom: 6px;">Set warnings</div>
						{#each ggResults.warnings as warning}
							<p style="display: flex; align-items: center; gap: 8px; font-size: 12.5px; color: var(--ink-soft); margin: 2px 0;">
								<span class="pw-dot pw-dot-warn"></span>{warning}
							</p>
						{/each}
					</div>
				{/if}
			</div>
		{/if}
	</div>
{/if}
