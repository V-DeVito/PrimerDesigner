<script>
	import '../app.css';

	let { children } = $props();
</script>

<div class="min-h-screen flex flex-col">
	<!-- Bauhaus color strip -->
	<div class="flex h-[6px]">
		<div class="flex-1 bg-[var(--color-red)]"></div>
		<div class="flex-1 bg-[var(--color-yellow)]"></div>
		<div class="flex-1 bg-[var(--color-blue)]"></div>
		<div class="flex-1 bg-[var(--color-green)]"></div>
	</div>

	<!-- Header -->
	<header class="border-b-[2.5px] border-[var(--color-ink)] px-6 py-3 flex items-center justify-between bg-[var(--color-surface)]">
		<div class="flex items-center gap-3">
			<!-- 2x2 geometric logo mark -->
			<div class="grid grid-cols-2 w-5 h-5 border border-[var(--color-ink)]">
				<div class="bg-[var(--color-red)]"></div>
				<div class="bg-[var(--color-yellow)]"></div>
				<div class="bg-[var(--color-blue)]"></div>
				<div class="bg-[var(--color-green)]"></div>
			</div>
			<span class="t-display text-[15px] tracking-[0.02em]">
				PrimerDesigner
			</span>
		</div>
		<nav class="flex items-center gap-1">
			<a href="/" class="px-3 py-1 text-[12px] font-bold uppercase tracking-[0.08em]
						border-[2px] border-[var(--color-ink)] bg-[var(--color-ink)] text-[var(--color-cream)]
						hover:bg-[var(--color-red)] transition-colors">
				Analyze
			</a>
			<a href="/golden-gate" class="px-3 py-1 text-[12px] font-bold uppercase tracking-[0.08em]
						border-[2px] border-[var(--color-ink)] bg-[var(--color-surface)]
						hover:bg-[var(--color-ink)] hover:text-[var(--color-cream)] transition-colors">
				Golden Gate
			</a>
			<a href="https://github.com/chenzymes/primerdesignr" target="_blank" rel="noopener"
			   class="ml-2 px-2 py-1 text-[12px] border-[2px] border-[var(--color-ink)]
					  hover:bg-[var(--color-ink)] hover:text-[var(--color-cream)] transition-colors">
				&#9679;
			</a>
		</nav>
	</header>

	<!-- Content -->
	<main class="flex-1 px-6 py-6 max-w-[1100px] mx-auto w-full">
		{@render children()}
	</main>

	<!-- Footer -->
	<div class="flex h-[4px]">
		<div class="flex-1 bg-[var(--color-red)]"></div>
		<div class="flex-1 bg-[var(--color-yellow)]"></div>
		<div class="flex-1 bg-[var(--color-blue)]"></div>
		<div class="flex-1 bg-[var(--color-green)]"></div>
	</div>
	<footer class="border-t-[2px] border-[var(--color-ink)] px-6 py-3 text-[11px] text-[var(--color-text-muted)]
				   flex items-center justify-between bg-[var(--color-surface)]">
		<span>SantaLucia 2004 &middot; Mathews/Turner NN parameters</span>
		<span>Built by <a href="https://vicenzodevito.com" class="font-semibold hover:text-[var(--color-text)]">Vicenzo DeVito</a></span>
	</footer>
</div>
