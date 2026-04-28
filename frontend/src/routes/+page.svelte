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
ATGCGTACGTAGCTAGCTACGATCGATCGTACGTAGCTAGCTAGCGATCGATCGTACGTAGCTAGCTAGCATCGATCGATGCTAGCTAGCTAGCATCGATCGTACGTAGCTAGCTAGCATCGATCGATCGATCGATCGTACGTAGCTAGCTAGCATCGATCGATCGTACGATCGATCGTAGCTAGCTAGCTAGCTAGCATCGATCGTACGATCGATCGATCGATCGATCGTAGCTAGCTAGCATCGATCGATCGTAGCTAGCATCGATCGATCGTAGCTAGCTAGCATCGATCGATCGTAGCTAGCTAGCTAGCATCGATCG`;

	let active = $state('design');

	let naConc = $state(50);
	let mgConc = $state(0);
	let dntpConc = $state(0);
	let dnaConc = $state(250);
	let showConditions = $state(false);

	let templateInput = $state(sampleTemplate);
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

		try {
			designResults = await designPrimers({
				template: templateInput,
				product_min: Number(productMin),
				product_max: Number(productMax),
				primer_count: Number(primerCount),
				target_start: optionalNumber(targetStart),
				target_length: optionalNumber(targetLength),
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

	function severityBg(value, floor = -12) {
		if (value === null || value === undefined) return '';
		const s = Math.abs(Math.max(floor, Math.min(0, value)) / floor);
		return `background: rgba(0,0,0,${(s * 0.12).toFixed(3)})`;
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

	function tabClass(name) {
		return active === name
			? 'border-gray-900 text-gray-900'
			: 'border-transparent text-gray-400 hover:text-gray-900';
	}
</script>

<svelte:head><title>PrimerDesigner</title></svelte:head>

<section class="mb-6">
	<div class="flex flex-col gap-3 sm:flex-row sm:items-end sm:justify-between">
		<div>
			<p class="text-[11px] uppercase tracking-[0.16em] text-gray-400">Primer workbench</p>
			<h1 class="mt-1 text-2xl font-semibold tracking-tight text-gray-900">Design, analyze, and explain primer sets.</h1>
		</div>
		<button
			class="self-start border border-gray-200 px-3 py-1.5 text-[12px] text-gray-500 transition-colors hover:border-gray-400 hover:text-gray-900"
			onclick={() => showConditions = !showConditions}
		>
			{showConditions ? 'Hide' : 'Show'} conditions
		</button>
	</div>

	{#if showConditions}
		<div class="mt-4 grid gap-3 border border-gray-200 bg-gray-50 p-3 text-[12px] text-gray-600 sm:grid-cols-4">
			<label>
				<span class="block text-[10px] uppercase tracking-wider text-gray-400">Na+ mM</span>
				<input class="mt-1 w-full border border-gray-200 bg-white px-2 py-1 font-mono" type="number" min="0" max="1000" step="10" bind:value={naConc} />
			</label>
			<label>
				<span class="block text-[10px] uppercase tracking-wider text-gray-400">Mg2+ mM</span>
				<input class="mt-1 w-full border border-gray-200 bg-white px-2 py-1 font-mono" type="number" min="0" max="100" step="0.5" bind:value={mgConc} />
			</label>
			<label>
				<span class="block text-[10px] uppercase tracking-wider text-gray-400">dNTP mM</span>
				<input class="mt-1 w-full border border-gray-200 bg-white px-2 py-1 font-mono" type="number" min="0" max="20" step="0.1" bind:value={dntpConc} />
			</label>
			<label>
				<span class="block text-[10px] uppercase tracking-wider text-gray-400">Primer nM</span>
				<input class="mt-1 w-full border border-gray-200 bg-white px-2 py-1 font-mono" type="number" min="1" max="10000" step="50" bind:value={dnaConc} />
			</label>
		</div>
	{/if}
</section>

<nav class="mb-6 flex border-b border-gray-200 text-sm font-medium">
	<button class={`border-b-2 px-4 py-2 transition-colors ${tabClass('design')}`} onclick={() => active = 'design'}>Design</button>
	<button class={`border-b-2 px-4 py-2 transition-colors ${tabClass('analyze')}`} onclick={() => active = 'analyze'}>Analyze</button>
	<button class={`border-b-2 px-4 py-2 transition-colors ${tabClass('golden')}`} onclick={() => active = 'golden'}>Golden Gate</button>
</nav>

{#if active === 'design'}
	<section class="grid gap-5 lg:grid-cols-[minmax(0,1fr)_320px]">
		<div>
			<div class="mb-2 flex items-center justify-between">
				<h2 class="text-sm font-semibold text-gray-900">Template sequence</h2>
				<span class="font-mono text-[11px] text-gray-400">{templateLength} bp</span>
			</div>
			<textarea
				class="h-60 w-full resize-y border border-gray-200 bg-gray-50 p-4 font-mono text-[12px] leading-6 focus:border-gray-400 focus:outline-none"
				bind:value={templateInput}
				spellcheck="false"
			></textarea>
		</div>

		<aside class="border border-gray-200 bg-white p-4">
			<h2 class="text-sm font-semibold text-gray-900">Design constraints</h2>
			<div class="mt-4 grid grid-cols-2 gap-3 text-[12px]">
				<label>
					<span class="block text-[10px] uppercase tracking-wider text-gray-400">Min product</span>
					<input class="mt-1 w-full border border-gray-200 px-2 py-1.5 font-mono" type="number" min="40" max="5000" bind:value={productMin} />
				</label>
				<label>
					<span class="block text-[10px] uppercase tracking-wider text-gray-400">Max product</span>
					<input class="mt-1 w-full border border-gray-200 px-2 py-1.5 font-mono" type="number" min="40" max="5000" bind:value={productMax} />
				</label>
				<label>
					<span class="block text-[10px] uppercase tracking-wider text-gray-400">Candidates</span>
					<input class="mt-1 w-full border border-gray-200 px-2 py-1.5 font-mono" type="number" min="1" max="20" bind:value={primerCount} />
				</label>
				<div></div>
				<label>
					<span class="block text-[10px] uppercase tracking-wider text-gray-400">Target start</span>
					<input class="mt-1 w-full border border-gray-200 px-2 py-1.5 font-mono" type="number" min="0" placeholder="optional" bind:value={targetStart} />
				</label>
				<label>
					<span class="block text-[10px] uppercase tracking-wider text-gray-400">Target length</span>
					<input class="mt-1 w-full border border-gray-200 px-2 py-1.5 font-mono" type="number" min="1" placeholder="optional" bind:value={targetLength} />
				</label>
			</div>
			<button
				class="mt-4 w-full bg-gray-900 px-4 py-2 text-sm font-medium text-white transition-colors hover:bg-black disabled:cursor-not-allowed disabled:opacity-30"
				disabled={designLoading || templateLength < 40}
				onclick={runDesign}
			>
				{designLoading ? 'Designing...' : 'Design primers'}
			</button>
			{#if designError}<p class="mt-3 text-[12px] text-red-600">{designError}</p>{/if}
		</aside>
	</section>

	{#if designResults}
		<section class="mt-6">
			<div class="mb-3 flex items-center justify-between">
				<div>
					<h2 class="text-sm font-semibold text-gray-900">Primer candidates</h2>
					<p class="text-[12px] text-gray-400">{designResults.candidates.length} returned for {designResults.template_length} bp template</p>
				</div>
				<button class="border border-gray-200 px-3 py-1.5 text-[12px] text-gray-500 hover:border-gray-400 hover:text-gray-900"
					onclick={() => downloadCsv('primerdesign-candidates.csv', exportDesignCSV(designResults))}>
					Export CSV
				</button>
			</div>

			{#if designResults.warnings.length || designResults.candidates.length === 0}
				<div class="mb-4 border border-gray-200 bg-gray-50 p-4">
					{#if designResults.warnings.length}
						{#each designResults.warnings as warning}
							<p class="text-[12px] text-gray-600">{warning}</p>
						{/each}
					{/if}
					{#if designResults.candidates.length === 0 && Object.keys(designResults.primer3_explain).length}
						<div class="mt-3 grid gap-2 text-[11px] text-gray-500 md:grid-cols-2">
							{#each Object.entries(designResults.primer3_explain) as [key, value]}
								<p><span class="font-mono text-gray-400">{key}</span>: {value}</p>
							{/each}
						</div>
					{/if}
				</div>
			{/if}

			{#if designResults.candidates.length}
			<div class="overflow-x-auto border border-gray-200">
				<table class="w-full border-collapse text-[13px]">
					<thead class="bg-gray-50 text-[10px] uppercase tracking-wider text-gray-400">
						<tr>
							<th class="px-3 py-2 text-left">Rank</th>
							<th class="px-3 py-2 text-left">Forward</th>
							<th class="px-3 py-2 text-left">Reverse</th>
							<th class="px-3 py-2 text-right">Product</th>
							<th class="px-3 py-2 text-right">Tm diff</th>
							<th class="px-3 py-2 text-right">Dimer</th>
						</tr>
					</thead>
					<tbody>
						{#each designResults.candidates as c}
							<tr class="cursor-pointer border-t border-gray-100 hover:bg-gray-50" onclick={() => expandedCandidate = expandedCandidate === c.rank ? null : c.rank}>
								<td class="px-3 py-2 font-mono">{c.rank}</td>
								<td class="px-3 py-2 font-mono text-[12px]">{c.forward.sequence}</td>
								<td class="px-3 py-2 font-mono text-[12px]">{c.reverse.sequence}</td>
								<td class="px-3 py-2 text-right font-mono">{c.product_size}</td>
								<td class="px-3 py-2 text-right font-mono">{c.tm_difference.toFixed(1)}</td>
								<td class="px-3 py-2 text-right font-mono" style={severityBg(c.heterodimer.dg)}>{c.heterodimer.dg.toFixed(2)}</td>
							</tr>
							{#if expandedCandidate === c.rank}
								<tr>
									<td colspan="6" class="border-t border-gray-100 bg-gray-50 px-4 py-4">
										<div class="grid gap-4 md:grid-cols-3">
											<div>
												<p class="text-[10px] uppercase tracking-wider text-gray-400">Forward coordinates</p>
												<p class="mt-1 font-mono text-[12px]">{c.forward_coords.start}-{c.forward_coords.end} ({c.forward_coords.strand})</p>
												<p class="mt-2 text-[12px] text-gray-500">Tm {c.forward.tm} C, GC {(c.forward.gc_content * 100).toFixed(0)}%</p>
											</div>
											<div>
												<p class="text-[10px] uppercase tracking-wider text-gray-400">Reverse coordinates</p>
												<p class="mt-1 font-mono text-[12px]">{c.reverse_coords.start}-{c.reverse_coords.end} ({c.reverse_coords.strand})</p>
												<p class="mt-2 text-[12px] text-gray-500">Tm {c.reverse.tm} C, GC {(c.reverse.gc_content * 100).toFixed(0)}%</p>
											</div>
											<div>
												<p class="text-[10px] uppercase tracking-wider text-gray-400">Actions</p>
												<button class="mt-1 border border-gray-200 px-3 py-1.5 text-[12px] hover:border-gray-400" onclick={() => copyPrimerPair(c)}>Copy pair</button>
											</div>
										</div>
										<div class="mt-4 grid gap-4 md:grid-cols-2">
											<div>
												<p class="text-[10px] uppercase tracking-wider text-gray-400">Explanation</p>
												{#each c.explanations as item}
													<p class="mt-1 text-[12px] text-gray-600">{item}</p>
												{/each}
											</div>
											<div>
												<p class="text-[10px] uppercase tracking-wider text-gray-400">Warnings</p>
												{#if c.warnings.length}
													{#each c.warnings as warning}
														<p class="mt-1 text-[12px] text-gray-600">{warning}</p>
													{/each}
												{:else}
													<p class="mt-1 text-[12px] text-gray-400">No warnings.</p>
												{/if}
											</div>
										</div>
									</td>
								</tr>
							{/if}
						{/each}
					</tbody>
				</table>
			</div>
			{/if}
		</section>
	{/if}
{/if}

{#if active === 'analyze'}
	<section>
		<div class="mb-2 flex items-center justify-between">
			<h2 class="text-sm font-semibold text-gray-900">Primer set</h2>
			<span class="font-mono text-[11px] text-gray-400">{analyzePrimerCount} detected</span>
		</div>
		<textarea
			class="h-44 w-full resize-y border border-gray-200 bg-gray-50 p-4 font-mono text-[12px] leading-6 focus:border-gray-400 focus:outline-none"
			bind:value={analyzeInput}
			spellcheck="false"
			placeholder={"F_primer\tATGCGATCG...\nR_primer\tCGATCGATC..."}
		></textarea>
		<div class="mt-3 flex items-center gap-3">
			<button class="bg-gray-900 px-4 py-2 text-sm font-medium text-white hover:bg-black disabled:cursor-not-allowed disabled:opacity-30"
				disabled={analyzeLoading || analyzePrimerCount === 0}
				onclick={runAnalyze}>
				{analyzeLoading ? 'Analyzing...' : 'Analyze set'}
			</button>
			{#if analyzeResults}
				<button class="border border-gray-200 px-3 py-2 text-[12px] text-gray-500 hover:border-gray-400 hover:text-gray-900"
					onclick={() => downloadCsv('primerdesign-analysis.csv', exportCSV(analyzeResults))}>
					Export CSV
				</button>
			{/if}
		</div>
		{#if analyzeError}<p class="mt-3 text-[12px] text-red-600">{analyzeError}</p>{/if}
	</section>

	{#if analyzeResults}
		<section class="mt-6">
			<div class="mb-4 grid grid-cols-3 gap-px border border-gray-200 bg-gray-200">
				<div class="bg-white px-4 py-3">
					<p class="text-[10px] uppercase tracking-wider text-gray-400">Warnings</p>
					<p class="mt-1 font-mono text-lg">{warningCount}<span class="text-sm text-gray-400">/{primerNames.length}</span></p>
				</div>
				<div class="bg-white px-4 py-3">
					<p class="text-[10px] uppercase tracking-wider text-gray-400">Tm spread</p>
					<p class="mt-1 font-mono text-lg">{analyzeResults.tm_spread.toFixed(1)} C</p>
				</div>
				<div class="bg-white px-4 py-3">
					<p class="text-[10px] uppercase tracking-wider text-gray-400">Range</p>
					<p class="mt-1 font-mono text-lg">{tmRange ? `${tmRange.min.toFixed(1)}-${tmRange.max.toFixed(1)} C` : '-'}</p>
				</div>
			</div>

			<div class="overflow-x-auto border border-gray-200">
				<table class="w-full border-collapse text-[13px]">
					<thead class="bg-gray-50 text-[10px] uppercase tracking-wider text-gray-400">
						<tr>
							<th class="px-3 py-2 text-left">Name</th>
							<th class="px-3 py-2 text-left">Sequence</th>
							<th class="px-3 py-2 text-right">Tm</th>
							<th class="px-3 py-2 text-right">GC</th>
							<th class="px-3 py-2 text-right">Hairpin</th>
							<th class="px-3 py-2 text-right">Self-dimer</th>
						</tr>
					</thead>
					<tbody>
						{#each Object.entries(analyzeResults.primers) as [name, p]}
							{@const isExpanded = expandedPrimer === name}
							<tr class="cursor-pointer border-t border-gray-100 hover:bg-gray-50" onclick={() => expandedPrimer = isExpanded ? null : name}>
								<td class="px-3 py-2 font-medium">{p.warnings.length ? '!' : ''} {name}</td>
								<td class="max-w-[420px] truncate px-3 py-2 font-mono text-[12px]" title={p.sequence}>{p.sequence}</td>
								<td class="px-3 py-2 text-right font-mono">{p.tm}</td>
								<td class="px-3 py-2 text-right font-mono">{(p.gc_content * 100).toFixed(0)}%</td>
								<td class="px-3 py-2 text-right font-mono" style={severityBg(p.hairpin.dg_santalucia, -6)}>{p.hairpin.dg_santalucia}</td>
								<td class="px-3 py-2 text-right font-mono" style={severityBg(p.homodimer.dg)}>{p.homodimer.dg}</td>
							</tr>
							{#if isExpanded}
								<tr>
									<td colspan="6" class="border-t border-gray-100 bg-gray-50 px-4 py-3">
										<div class="grid gap-4 text-[12px] md:grid-cols-3">
											<p><span class="text-gray-400">Mathews hairpin</span> <span class="font-mono">{p.hairpin.dg_mathews}</span></p>
											<p><span class="text-gray-400">Dot bracket</span> <span class="font-mono">{p.hairpin.dot_bracket || '-'}</span></p>
											<p><span class="text-gray-400">Engines disagree</span> <span class="font-mono">{p.hairpin.engines_disagree ? 'yes' : 'no'}</span></p>
										</div>
										{#if p.warnings.length}
											<div class="mt-3">
												{#each p.warnings as warning}
													<p class="text-[12px] text-gray-600">{warning}</p>
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

			<div class="mt-6 overflow-x-auto">
				<p class="mb-2 text-[10px] uppercase tracking-wider text-gray-400">Cross-dimer matrix, delta G kcal/mol</p>
				<table class="w-full table-fixed border-collapse font-mono text-[12px]">
					<thead>
						<tr>
							<th class="w-24"></th>
							{#each primerNames as name}
								<th class="truncate pb-1 text-center text-[10px] font-normal uppercase tracking-wider text-gray-400">{name}</th>
							{/each}
						</tr>
					</thead>
					<tbody>
						{#each primerNames as row, i}
							<tr>
								<td class="truncate pr-2 text-right text-[10px] uppercase tracking-wider text-gray-400">{row}</td>
								{#each primerNames as col, j}
									{#if i === j}
										{@const selfDg = analyzeResults.primers[row].homodimer.dg}
										<td class="border border-gray-200 py-1.5 text-center" style={severityBg(selfDg)}>{selfDg}</td>
									{:else}
										{@const cd = getCrossDimer(row, col)}
										<td class="border border-gray-200 py-1.5 text-center" style={cd ? severityBg(cd.dg) : ''}>{cd ? cd.dg.toFixed(1) : ''}</td>
									{/if}
								{/each}
							</tr>
						{/each}
					</tbody>
				</table>
			</div>
		</section>
	{/if}
{/if}

{#if active === 'golden'}
	<section class="grid gap-5 lg:grid-cols-[minmax(0,1fr)_320px]">
		<div>
			<h2 class="mb-2 text-sm font-semibold text-gray-900">Overhangs</h2>
			<textarea class="h-36 w-full resize-y border border-gray-200 bg-gray-50 p-4 font-mono text-[12px] leading-6 focus:border-gray-400 focus:outline-none"
				bind:value={overhangInput}
				spellcheck="false"></textarea>
			<h2 class="mb-2 mt-4 text-sm font-semibold text-gray-900">Sequences to scan for internal sites</h2>
			<textarea class="h-36 w-full resize-y border border-gray-200 bg-gray-50 p-4 font-mono text-[12px] leading-6 focus:border-gray-400 focus:outline-none"
				bind:value={ggSequences}
				placeholder={">insert\nGGTCTC..."}
				spellcheck="false"></textarea>
		</div>
		<aside class="border border-gray-200 bg-white p-4">
			<label class="block text-[12px]">
				<span class="block text-[10px] uppercase tracking-wider text-gray-400">Type IIS enzyme</span>
				<select class="mt-1 w-full border border-gray-200 bg-white px-2 py-2" bind:value={enzyme}>
					<option>BsaI</option>
					<option>BsmBI</option>
					<option>BbsI</option>
					<option>SapI</option>
					<option>BpiI</option>
				</select>
			</label>
			<button class="mt-4 w-full bg-gray-900 px-4 py-2 text-sm font-medium text-white hover:bg-black disabled:opacity-30"
				disabled={ggLoading || parseOverhangs(overhangInput).length < 2}
				onclick={runGoldenGate}>
				{ggLoading ? 'Checking...' : 'Check overhangs'}
			</button>
			{#if ggError}<p class="mt-3 text-[12px] text-red-600">{ggError}</p>{/if}
		</aside>
	</section>

	{#if ggResults}
		<section class="mt-6">
			<h2 class="mb-3 text-sm font-semibold text-gray-900">Golden Gate report</h2>
			<div class="overflow-x-auto border border-gray-200">
				<table class="w-full border-collapse text-[13px]">
					<thead class="bg-gray-50 text-[10px] uppercase tracking-wider text-gray-400">
						<tr>
							<th class="px-3 py-2 text-left">Overhang</th>
							<th class="px-3 py-2 text-left">Reverse complement</th>
							<th class="px-3 py-2 text-right">GC</th>
							<th class="px-3 py-2 text-left">Fidelity set</th>
							<th class="px-3 py-2 text-left">Warnings</th>
						</tr>
					</thead>
					<tbody>
						{#each ggResults.overhangs as oh}
							<tr class="border-t border-gray-100">
								<td class="px-3 py-2 font-mono">{oh.overhang}</td>
								<td class="px-3 py-2 font-mono">{oh.reverse_complement}</td>
								<td class="px-3 py-2 text-right font-mono">{(oh.gc_content * 100).toFixed(0)}%</td>
								<td class="px-3 py-2">{oh.is_in_fidelity_set ? 'NEB high fidelity' : 'custom'}</td>
								<td class="px-3 py-2 text-gray-600">{oh.warnings.join('; ') || '-'}</td>
							</tr>
						{/each}
					</tbody>
				</table>
			</div>
			{#if ggResults.warnings.length}
				<div class="mt-4 border border-gray-200 bg-gray-50 p-4">
					<p class="text-[10px] uppercase tracking-wider text-gray-400">Warnings</p>
					{#each ggResults.warnings as warning}
						<p class="mt-1 text-[12px] text-gray-600">{warning}</p>
					{/each}
				</div>
			{/if}
		</section>
	{/if}
{/if}
