<script>
	import '../app.css';

	let { children } = $props();
</script>

<div class="min-h-screen flex flex-col">
	<header class="border-b border-gray-200 px-6 py-3 flex items-center justify-between">
		<a href="/" class="text-sm font-semibold tracking-tight text-gray-900">PrimerDesigner</a>
		<nav class="flex items-center gap-4 text-[13px]">
			<a href="/" class="text-gray-500 hover:text-gray-900 transition-colors">Analyze</a>
			<a href="/golden-gate" class="text-gray-500 hover:text-gray-900 transition-colors">Golden Gate</a>
			<a href="https://github.com/V-DeVito/PrimerDesigner" target="_blank" rel="noopener"
			   class="text-gray-400 hover:text-gray-900 transition-colors">
				GitHub
			</a>
		</nav>
	</header>

	<main class="flex-1 max-w-[960px] w-full mx-auto px-6 py-8">
		{@render children()}
	</main>

	<footer class="border-t border-gray-200 px-6 py-3 text-[11px] text-gray-400 flex items-center justify-between">
		<span>SantaLucia 2004 · Mathews/Turner NN parameters</span>
		<span>Built by <a href="https://vicenzodevito.com" class="hover:text-gray-900">Vicenzo DeVito</a></span>
	</footer>
</div>
