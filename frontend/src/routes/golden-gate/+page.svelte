<script>
	import { validateGoldenGate } from '$lib/api.js';

	// ── State ───────────────────────────────────────────
	let overhangInput = $state('AATG, TTCG, GCAA, CCTA');
	let enzyme = $state('BsaI');
	let seqInput = $state('');
	let primerInput = $state('');
	let loading = $state(false);
	let error = $state('');
	let results = $state(null);
	let showAdvanced = $state(false);
	let toast = $state('');
	let toastKey = $state(0);

	const enzymes = ['BsaI', 'BbsI', 'BsmBI', 'SapI', 'BtgZI'];

	// ── Derived ─────────────────────────────────────────
	let parsedOverhangs = $derived.by(() => {
		return overhangInput
			.split(/[\s,]+/)
			.map(s => s.trim().toUpperCase())
			.filter(s => /^[ATGC]{4}$/.test(s));
	});

	let overhangCount = $derived(parsedOverhangs.length);

	// ── Actions ─────────────────────────────────────────
	async function validate() {
		if (parsedOverhangs.length === 0) {
			error = 'Enter at least one 4bp overhang.';
			return;
		}
		loading = true;
		error = '';
		results = null;

		// Parse sequences if provided
		let sequences = undefined;
		if (seqInput.trim()) {
			sequences = {};
			const lines = seqInput.trim().split('\n');
			for (let i = 0; i < lines.length; i++) {
				const line = lines[i].trim();
				if (!line) continue;
				if (line.startsWith('>')) {
					const name = line.slice(1).trim();
					i++;
					if (i < lines.length) sequences[name] = lines[i].trim().toUpperCase();
				} else {
					const match = line.match(/^(.+?)[\t, ]+([ATGCatgc]+)$/);
					if (match) sequences[match[1].trim()] = match[2].toUpperCase();
				}
			}
			if (Object.keys(sequences).length === 0) sequences = undefined;
		}

		// Parse primers if provided
		let primers = undefined;
		if (primerInput.trim()) {
			primers = {};
			const lines = primerInput.trim().split('\n');
			for (const line of lines) {
				const match = line.trim().match(/^(.+?)[\t, ]+([ATGCatgc]+)$/);
				if (match) primers[match[1].trim()] = match[2].toUpperCase();
			}
			if (Object.keys(primers).length === 0) primers = undefined;
		}

		try {
			results = await validateGoldenGate(parsedOverhangs, enzyme, sequences, primers);
		} catch (e) {
			error = e.message;
		} finally {
			loading = false;
		}
	}

	function showToast(msg) {
		toast = msg;
		toastKey++;
		setTimeout(() => { toast = ''; }, 2200);
	}

	function gcColor(gc) {
		if (gc >= 0.25 && gc <= 0.75) return 'val-ok';
		return 'val-warn';
	}

	function handleKeydown(e) {
		if ((e.metaKey || e.ctrlKey) && e.key === 'Enter') {
			e.preventDefault();
			validate();
		}
	}
</script>

<svelte:head>
	<title>Golden Gate — primerdesignr</title>
</svelte:head>

<svelte:window on:keydown={handleKeydown} />

<div class="max-w-4xl mx-auto px-5 py-6">
	<!-- Header -->
	<div class="mb-6">
		<h1 class="text-xl font-semibold tracking-tight">Golden Gate Assembly</h1>
		<p class="text-sm text-[var(--color-text-secondary)] mt-1">
			Validate overhang sets for uniqueness, GC content, palindromes, and internal enzyme sites.
		</p>
	</div>

	<!-- Input Section -->
	<section class="mb-8">
		<div class="card !p-5">
			<!-- Overhangs -->
			<div class="mb-4">
				<label class="text-xs font-semibold text-[var(--color-text-secondary)] uppercase tracking-wider block mb-1.5">
					Overhangs (4bp each)
				</label>
				<input
					type="text"
					class="input input-mono w-full"
					bind:value={overhangInput}
					placeholder="AATG, TTCG, GCAA, CCTA"
					spellcheck="false"
				/>
				<div class="text-[10px] text-[var(--color-text-muted)] mt-1.5">
					Comma or space-separated. &thinsp;
					{#if overhangCount > 0}
						<span class="font-mono">{overhangCount}</span> overhang{overhangCount !== 1 ? 's' : ''} detected.
					{/if}
				</div>
			</div>

			<!-- Enzyme selector -->
			<div class="mb-4">
				<label class="text-xs font-semibold text-[var(--color-text-secondary)] uppercase tracking-wider block mb-1.5">
					Enzyme
				</label>
				<div class="flex gap-1.5">
					{#each enzymes as enz}
						<button
							class="px-3 py-1.5 rounded-md text-sm font-mono font-medium transition-all
							       {enzyme === enz
							         ? 'bg-[var(--color-primary)] text-white'
							         : 'bg-[var(--color-surface-raised)] text-[var(--color-text-secondary)] hover:text-[var(--color-text)] hover:bg-[var(--color-surface-overlay)]'}"
							onclick={() => enzyme = enz}
						>
							{enz}
						</button>
					{/each}
				</div>
			</div>

			<!-- Advanced: sequences + primers -->
			<button
				class="text-xs text-[var(--color-text-muted)] hover:text-[var(--color-text-secondary)]
				       flex items-center gap-1 mb-3 transition-colors"
				onclick={() => showAdvanced = !showAdvanced}
			>
				<svg class="w-3 h-3 transition-transform {showAdvanced ? 'rotate-180' : ''}"
				     fill="none" stroke="currentColor" viewBox="0 0 24 24">
					<path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M19 9l-7 7-7-7"/>
				</svg>
				Scan for internal sites (optional)
			</button>

			{#if showAdvanced}
				<div class="grid gap-4 sm:grid-cols-2 mb-4 fade-in">
					<div>
						<label class="text-[10px] font-medium text-[var(--color-text-muted)] uppercase tracking-wider block mb-1">
							Part sequences
						</label>
						<textarea
							class="seq-input text-xs"
							rows="3"
							bind:value={seqInput}
							placeholder={"Part1\tATGCGATCG...\nPart2\tCGATCGATC..."}
							spellcheck="false"
						></textarea>
					</div>
					<div>
						<label class="text-[10px] font-medium text-[var(--color-text-muted)] uppercase tracking-wider block mb-1">
							Primers to scan
						</label>
						<textarea
							class="seq-input text-xs"
							rows="3"
							bind:value={primerInput}
							placeholder={"F_primer\tATGCGATCG...\nR_primer\tCGATCGATC..."}
							spellcheck="false"
						></textarea>
					</div>
				</div>
			{/if}

			<!-- Validate button -->
			<div class="flex justify-end gap-2 items-center">
				<span class="text-[10px] text-[var(--color-text-muted)]">
					{#if !loading}Cmd+Enter{/if}
				</span>
				<button
					onclick={validate}
					disabled={loading || overhangCount === 0}
					class="btn btn-primary"
				>
					{#if loading}
						<span class="spinner"></span>
						Validating...
					{:else}
						Validate
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

	<!-- Results -->
	{#if results}
		<!-- Summary -->
		<section class="mb-6 fade-in">
			<div class="grid grid-cols-2 gap-3">
				<div class="card flex items-center gap-3">
					{#if results.all_unique}
						<div class="w-8 h-8 rounded-lg bg-[var(--color-ok-soft)] flex items-center justify-center">
							<svg class="w-4 h-4 text-[var(--color-ok)]" fill="none" stroke="currentColor" viewBox="0 0 24 24">
								<path stroke-linecap="round" stroke-linejoin="round" stroke-width="2.5" d="M5 13l4 4L19 7"/>
							</svg>
						</div>
						<div>
							<div class="text-sm font-semibold text-[var(--color-ok)]">All Unique</div>
							<div class="text-[11px] text-[var(--color-text-muted)]">No collisions detected</div>
						</div>
					{:else}
						<div class="w-8 h-8 rounded-lg bg-[var(--color-danger-soft)] flex items-center justify-center">
							<svg class="w-4 h-4 text-[var(--color-danger)]" fill="none" stroke="currentColor" viewBox="0 0 24 24">
								<path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M6 18L18 6M6 6l12 12"/>
							</svg>
						</div>
						<div>
							<div class="text-sm font-semibold text-[var(--color-danger)]">Collisions</div>
							<div class="text-[11px] text-[var(--color-text-muted)]">Non-unique overhangs found</div>
						</div>
					{/if}
				</div>

				<div class="card">
					<div class="text-[11px] text-[var(--color-text-muted)] uppercase tracking-wider font-medium mb-1">Enzyme</div>
					<div class="text-lg font-semibold font-mono">{results.enzyme}</div>
					{#if results.internal_sites?.length > 0}
						<div class="text-[10px] text-[var(--color-danger)] font-medium mt-0.5">
							{results.internal_sites.length} internal site{results.internal_sites.length > 1 ? 's' : ''} found
						</div>
					{/if}
				</div>
			</div>
		</section>

		<!-- Overhang cards -->
		<section class="mb-6">
			<h2 class="text-xs font-semibold text-[var(--color-text-muted)] uppercase tracking-wider mb-3">
				Overhangs
			</h2>
			<div class="grid gap-2 sm:grid-cols-2">
				{#each results.overhangs as oh, idx}
					{@const hasIssues = oh.warnings.length > 0 || oh.is_palindromic}
					<div class="card {hasIssues ? 'card-warn' : 'card-ok'} fade-in"
					     style="animation-delay: {idx * 40}ms">
						<div class="flex items-center justify-between mb-2">
							<span class="font-mono text-base font-semibold tracking-widest">
								{oh.overhang}
							</span>
							<div class="flex gap-1.5">
								{#if oh.is_palindromic}
									<span class="badge badge-warn">palindrome</span>
								{/if}
								{#if oh.is_in_fidelity_set}
									<span class="badge badge-ok">NEB fidelity</span>
								{:else}
									<span class="badge badge-warn">non-standard</span>
								{/if}
							</div>
						</div>
						<div class="flex gap-6 text-xs">
							<div class="metric">
								<span class="metric-label">RC</span>
								<span class="metric-value font-mono text-[var(--color-text-secondary)]">{oh.reverse_complement}</span>
							</div>
							<div class="metric">
								<span class="metric-label">GC</span>
								<span class="metric-value {gcColor(oh.gc_content)}">{(oh.gc_content * 100).toFixed(0)}%</span>
							</div>
						</div>
						{#if oh.warnings.length > 0}
							<div class="mt-2 pt-2 border-t border-[var(--color-border-subtle)]">
								{#each oh.warnings as w}
									<div class="text-xs text-[var(--color-warn)] flex items-start gap-1.5">
										<svg class="w-3 h-3 shrink-0 mt-0.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
											<path stroke-linecap="round" stroke-linejoin="round" stroke-width="2"
											      d="M12 9v2m0 4h.01"/>
										</svg>
										{w}
									</div>
								{/each}
							</div>
						{/if}
					</div>
				{/each}
			</div>
		</section>

		<!-- Internal sites -->
		{#if results.internal_sites?.length > 0}
			<section class="mb-6 fade-in">
				<h2 class="text-xs font-semibold text-[var(--color-text-muted)] uppercase tracking-wider mb-3">
					Internal Enzyme Sites
				</h2>
				<div class="card card-danger">
					{#each results.internal_sites as site}
						<div class="flex items-center gap-3 text-sm py-1">
							<span class="font-mono text-[var(--color-danger)]">pos {site.position}</span>
							<span class="text-[var(--color-text-secondary)]">in</span>
							<span class="font-medium">{site.sequence}</span>
						</div>
					{/each}
				</div>
			</section>
		{/if}

		<!-- Global warnings -->
		{#if results.warnings?.length > 0}
			<section class="mb-6 fade-in">
				<div class="card card-warn">
					{#each results.warnings as w}
						<div class="flex items-start gap-2 text-sm text-[var(--color-warn)] py-0.5">
							<svg class="w-4 h-4 shrink-0 mt-0.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
								<path stroke-linecap="round" stroke-linejoin="round" stroke-width="2"
								      d="M12 9v2m0 4h.01m-6.938 4h13.856c1.54 0 2.502-1.667 1.732-2.5L13.732 4.5c-.77-.833-2.694-.833-3.464 0L3.34 16.5c-.77.833.192 2.5 1.732 2.5z"/>
							</svg>
							{w}
						</div>
					{/each}
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
