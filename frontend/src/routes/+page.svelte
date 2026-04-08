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

	// Severity: 0 (benign) → 1 (bad). Returns opacity for black fill.
	function dgSeverity(dg, floor = -12) {
		if (dg === null || dg === undefined) return 0;
		return Math.abs(Math.max(floor, Math.min(0, dg)) / floor);
	}

	function dgBg(dg, floor = -12) {
		if (dg === null || dg === undefined) return '';
		const s = dgSeverity(dg, floor);
		// Black at varying opacity — more severe = darker
		return `background: rgba(0,0,0,${(s * 0.12).toFixed(3)})`;
	}

	function handleKeydown(e) {
		if ((e.metaKey || e.ctrlKey) && e.key === 'Enter') { e.preventDefault(); analyze(); }
	}
</script>

<svelte:head><title>Analyze — PrimerDesigner</title></svelte:head>
<svelte:window on:keydown={handleKeydown} />

<!-- Input -->
<div class="mb-6">
	<textarea
		class="w-full font-mono text-[13px] leading-[1.8] p-4 border border-gray-200
			   bg-gray-50 resize-none focus:outline-none focus:border-gray-400 transition-colors"
		rows="6"
		bind:value={input}
		placeholder={"F_primer\tATGCGATCG...\nR_primer\tCGATCGATC..."}
		spellcheck="false"
	></textarea>
	<div class="flex items-center gap-3 mt-2">
		<button onclick={analyze} disabled={loading || primerCount === 0}
				class="px-4 py-1.5 text-[12px] font-medium bg-gray-900 text-white
					   hover:bg-black disabled:opacity-30 disabled:cursor-not-allowed transition-colors">
			{#if loading}Analyzing...{:else}Run{/if}
		</button>
		<button class="text-[12px] text-gray-400 hover:text-gray-700 transition-colors"
				onclick={() => showConditions = !showConditions}>
			{showConditions ? '−' : '+'} Conditions
		</button>
	</div>
	{#if showConditions}
		<div class="flex items-center gap-4 mt-2 text-[11px] text-gray-500">
			<label class="flex items-center gap-1">
				Na+
				<input type="number" bind:value={naConc} min="0" max="1000" step="10"
					class="w-12 text-center py-0.5 border border-gray-200 bg-white font-mono text-[11px]
						   focus:outline-none focus:border-gray-400" />
				mM
			</label>
			<label class="flex items-center gap-1">
				Mg2+
				<input type="number" bind:value={mgConc} min="0" max="100" step="0.5"
					class="w-12 text-center py-0.5 border border-gray-200 bg-white font-mono text-[11px]
						   focus:outline-none focus:border-gray-400" />
				mM
			</label>
			<label class="flex items-center gap-1">
				DNA
				<input type="number" bind:value={dnaConc} min="1" max="10000" step="50"
					class="w-14 text-center py-0.5 border border-gray-200 bg-white font-mono text-[11px]
						   focus:outline-none focus:border-gray-400" />
				nM
			</label>
		</div>
	{/if}
	{#if error}
		<p class="mt-2 text-[12px] text-red-600">{error}</p>
	{/if}
</div>

{#if loading}
	<div class="py-20 text-center text-[13px] text-gray-400">Analyzing...</div>
{/if}

{#if results}
	<!-- Summary row -->
	<div class="grid grid-cols-3 gap-px bg-gray-200 border border-gray-200 mb-6">
		<div class="bg-white px-4 py-3">
			<div class="text-[10px] uppercase tracking-wider text-gray-400 mb-1">Warnings</div>
			<div class="text-lg font-semibold font-mono">
				{#if warningCount === 0}
					<span class="text-gray-300">None</span>
				{:else}
					{warningCount}<span class="text-gray-400 text-sm">/{primerNames.length}</span>
				{/if}
			</div>
		</div>
		<div class="bg-white px-4 py-3">
			<div class="text-[10px] uppercase tracking-wider text-gray-400 mb-1">Tm spread</div>
			<div class="text-lg font-semibold font-mono">{results.tm_spread.toFixed(1)}&deg;</div>
			{#if tmRange}
				<div class="text-[11px] text-gray-400 font-mono">{tmRange.min.toFixed(1)}&ndash;{tmRange.max.toFixed(1)}&deg;C</div>
			{/if}
		</div>
		<div class="bg-white px-4 py-3">
			<div class="text-[10px] uppercase tracking-wider text-gray-400 mb-1">GC range</div>
			{#if gcRange}
				<div class="text-lg font-semibold font-mono">{gcRange.min.toFixed(0)}&ndash;{gcRange.max.toFixed(0)}%</div>
				<div class="text-[11px] text-gray-400">target 40&ndash;60%</div>
			{/if}
		</div>
	</div>

	<!-- Primer table -->
	<table class="w-full text-[13px] border-collapse mb-6">
		<thead>
			<tr class="border-b border-gray-900 text-left">
				<th class="py-2 pr-4 text-[10px] uppercase tracking-wider text-gray-500 font-medium">Name</th>
				<th class="py-2 pr-4 text-[10px] uppercase tracking-wider text-gray-500 font-medium">Sequence</th>
				<th class="py-2 px-2 text-[10px] uppercase tracking-wider text-gray-500 font-medium text-right">Tm</th>
				<th class="py-2 px-2 text-[10px] uppercase tracking-wider text-gray-500 font-medium text-right">GC%</th>
				<th class="py-2 px-2 text-[10px] uppercase tracking-wider text-gray-500 font-medium text-right">Hairpin</th>
			</tr>
		</thead>
		<tbody>
			{#each Object.entries(results.primers) as [name, p], idx}
				{@const hasWarnings = p.warnings.length > 0}
				{@const isExpanded = expandedPrimer === name}
				<tr class="border-b border-gray-100 hover:bg-gray-50 cursor-pointer transition-colors"
					onclick={() => expandedPrimer = isExpanded ? null : name}
					onkeydown={(e) => e.key === 'Enter' && (expandedPrimer = isExpanded ? null : name)}
					role="button" tabindex="0">
					<td class="py-2.5 pr-4 font-medium whitespace-nowrap">
						{#if hasWarnings}<span class="text-gray-400 mr-1" title="{p.warnings.length} warning{p.warnings.length > 1 ? 's' : ''}">!</span>{/if}{name}
					</td>
					<td class="py-2.5 pr-4">
						<span class="seq whitespace-nowrap overflow-hidden text-ellipsis block max-w-[400px]" title={p.sequence}>
							{@html colorSeq(p.sequence)}
						</span>
					</td>
					<td class="py-2.5 px-2 font-mono text-right" style="{dgBg(-(Math.abs(p.tm - 60)), -10)}">{p.tm}</td>
					<td class="py-2.5 px-2 font-mono text-right" style="{dgBg(-(Math.abs(p.gc_content * 100 - 50)), -20)}">{(p.gc_content * 100).toFixed(0)}</td>
					<td class="py-2.5 px-2 font-mono text-right" style="{dgBg(p.hairpin.dg_santalucia, -6)}">{p.hairpin.dg_santalucia}</td>
				</tr>
				{#if isExpanded}
					<tr>
						<td colspan="5" class="py-3 px-4 bg-gray-50 border-b border-gray-100">
							<div class="grid grid-cols-4 gap-6 text-[12px] mb-2">
								<div>
									<span class="text-gray-400">Hairpin (SL)</span>
									<span class="font-mono ml-2">{p.hairpin.dg_santalucia}</span>
								</div>
								<div>
									<span class="text-gray-400">Hairpin (M)</span>
									<span class="font-mono ml-2">{p.hairpin.dg_mathews}</span>
									{#if p.hairpin.engines_disagree}<span class="text-gray-400 ml-1">*</span>{/if}
								</div>
								<div>
									<span class="text-gray-400">Self-dimer</span>
									<span class="font-mono ml-2">{p.homodimer.dg}</span>
								</div>
								<div>
									<span class="text-gray-400">3' dimer</span>
									<span class="font-mono ml-2">{p.homodimer.dg_3prime ?? '—'}</span>
								</div>
							</div>
							{#if hasWarnings}
								{#each p.warnings as w}
									<p class="text-[11px] text-gray-500 mt-1">! {w}</p>
								{/each}
							{/if}
						</td>
					</tr>
				{/if}
			{/each}
		</tbody>
	</table>

	<!-- Dimer matrix -->
	<div class="mb-6">
		<div class="flex items-baseline justify-between mb-2">
			<span class="text-[10px] uppercase tracking-wider text-gray-400 font-medium">Dimer matrix</span>
			<span class="text-[10px] text-gray-400">&#916;G kcal/mol</span>
		</div>
		<table class="w-full border-collapse text-[12px] font-mono" style="table-layout: fixed">
			<thead>
				<tr>
					<th class="w-[60px]"></th>
					{#each primerNames as name}
						<th class="text-[10px] text-gray-400 font-normal text-center pb-1.5 uppercase tracking-wider">{name}</th>
					{/each}
				</tr>
			</thead>
			<tbody>
				{#each primerNames as row, i}
					<tr>
						<td class="text-[10px] text-gray-400 text-right pr-2 uppercase tracking-wider">{row}</td>
						{#each primerNames as col, j}
							{#if i === j}
								{@const selfDg = results.primers[row].homodimer.dg}
								<td class="border border-gray-200 text-center py-1.5 font-medium"
									style="{dgBg(selfDg)}"
									title="{row} self-dimer: {selfDg}">
									{selfDg}
								</td>
							{:else}
								{@const cd = getCrossDimer(row, col)}
								<td class="border border-gray-200 text-center py-1.5"
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
{/if}
