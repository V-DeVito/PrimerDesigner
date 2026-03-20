<script>
	import '../app.css';
	import { page } from '$app/stores';

	let currentPath = $derived($page.url.pathname);
</script>

<div class="min-h-screen flex flex-col">
	<!-- Header -->
	<header class="sticky top-0 z-50 border-b border-[var(--color-border-subtle)]
	               bg-[var(--color-bg)]/80 backdrop-blur-xl">
		<div class="max-w-7xl mx-auto px-5 h-12 flex items-center justify-between">
			<!-- Logo -->
			<a href="/" class="flex items-center gap-2.5 group">
				<div class="w-6 h-6 rounded-md bg-[var(--color-primary)] flex items-center justify-center
				            text-white text-xs font-bold font-mono">P</div>
				<span class="text-sm font-semibold tracking-tight text-[var(--color-text)]
				             group-hover:text-[var(--color-primary)] transition-colors">
					primerdesignr
				</span>
				<span class="text-[10px] font-mono text-[var(--color-text-muted)] bg-[var(--color-surface-raised)]
				             px-1.5 py-0.5 rounded">v0.3</span>
			</a>

			<!-- Navigation tabs -->
			<nav class="flex items-center gap-1">
				<a href="/"
				   class="px-3 py-1.5 rounded-md text-sm font-medium transition-all
				          {currentPath === '/'
				            ? 'text-[var(--color-text)] bg-[var(--color-surface-raised)]'
				            : 'text-[var(--color-text-secondary)] hover:text-[var(--color-text)] hover:bg-[var(--color-surface)]'}">
					Analyze
				</a>
				<a href="/golden-gate"
				   class="px-3 py-1.5 rounded-md text-sm font-medium transition-all
				          {currentPath === '/golden-gate'
				            ? 'text-[var(--color-text)] bg-[var(--color-surface-raised)]'
				            : 'text-[var(--color-text-secondary)] hover:text-[var(--color-text)] hover:bg-[var(--color-surface)]'}">
					Golden Gate
				</a>
				<div class="w-px h-4 bg-[var(--color-border)] mx-2"></div>
				<a href="https://github.com/V-DeVito/PrimerDesigner" target="_blank" rel="noopener"
				   class="px-2 py-1.5 text-[var(--color-text-muted)] hover:text-[var(--color-text-secondary)] transition-colors">
					<svg class="w-4 h-4" fill="currentColor" viewBox="0 0 24 24">
						<path d="M12 0c-6.626 0-12 5.373-12 12 0 5.302 3.438 9.8 8.207 11.387.599.111.793-.261.793-.577v-2.234c-3.338.726-4.033-1.416-4.033-1.416-.546-1.387-1.333-1.756-1.333-1.756-1.089-.745.083-.729.083-.729 1.205.084 1.839 1.237 1.839 1.237 1.07 1.834 2.807 1.304 3.492.997.107-.775.418-1.305.762-1.604-2.665-.305-5.467-1.334-5.467-5.931 0-1.311.469-2.381 1.236-3.221-.124-.303-.535-1.524.117-3.176 0 0 1.008-.322 3.301 1.23.957-.266 1.983-.399 3.003-.404 1.02.005 2.047.138 3.006.404 2.291-1.552 3.297-1.23 3.297-1.23.653 1.653.242 2.874.118 3.176.77.84 1.235 1.911 1.235 3.221 0 4.609-2.807 5.624-5.479 5.921.43.372.823 1.102.823 2.222v3.293c0 .319.192.694.801.576 4.765-1.589 8.199-6.086 8.199-11.386 0-6.627-5.373-12-12-12z"/>
					</svg>
				</a>
			</nav>
		</div>
	</header>

	<!-- Content -->
	<main class="flex-1">
		<slot />
	</main>

	<!-- Footer -->
	<footer class="border-t border-[var(--color-border-subtle)] py-4">
		<div class="max-w-7xl mx-auto px-5 flex items-center justify-between text-[11px] text-[var(--color-text-muted)]">
			<span>SantaLucia 2004 + Mathews/Turner nearest-neighbor parameters</span>
			<span>
				Built by
				<a href="https://vicenzodevito.com" class="text-[var(--color-text-secondary)] hover:text-[var(--color-text)] transition-colors">
					Vicenzo DeVito
				</a>
			</span>
		</div>
	</footer>
</div>
