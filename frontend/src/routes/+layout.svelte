<script>
	import '../app.css';
</script>

<div class="min-h-screen flex flex-col">
	<!-- Header -->
	<header class="border-b border-[var(--color-surface-overlay)] px-6 py-3 flex items-center justify-between">
		<div class="flex items-center gap-3">
			<span class="font-mono text-sm font-semibold tracking-wider text-[var(--color-info)]">
				primerdesignr
			</span>
			<span class="text-xs text-[var(--color-text-dim)]">v0.1.0</span>
		</div>
		<nav class="flex gap-4 text-sm text-[var(--color-text-muted)]">
			<a href="/" class="hover:text-[var(--color-text)] transition-colors">Analyze</a>
			<a href="/golden-gate" class="hover:text-[var(--color-text)] transition-colors">Golden Gate</a>
			<a href="https://github.com/chenzymes/primerdesignr" target="_blank" rel="noopener"
			   class="hover:text-[var(--color-text)] transition-colors">GitHub</a>
		</nav>
	</header>

	<!-- Content -->
	<main class="flex-1 px-6 py-6 max-w-6xl mx-auto w-full">
		<slot />
	</main>

	<!-- Footer -->
	<footer class="border-t border-[var(--color-surface-overlay)] px-6 py-3 text-xs text-[var(--color-text-dim)] text-center">
		SantaLucia 2004 + Mathews/Turner · primer3-py · seqfold · Built by
		<a href="https://vicenzodevito.com" class="text-[var(--color-text-muted)] hover:text-[var(--color-text)]">Vicenzo DeVito</a>
	</footer>
</div>
