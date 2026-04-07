<script>
	import '../app.css';

	let { children } = $props();
</script>

<div class="bg-scene"></div>

<div class="min-h-screen flex flex-col relative">
	<header class="px-6 py-4 flex items-center justify-between">
		<a href="/" class="text-[14px] font-semibold tracking-tight text-white/90">PrimerDesigner</a>
		<nav class="flex items-center gap-4 text-[13px]">
			<a href="https://github.com/V-DeVito/PrimerDesigner" target="_blank" rel="noopener"
			   class="text-white/30 hover:text-white/70 transition-colors">
				GitHub
			</a>
		</nav>
	</header>

	<main class="flex-1 max-w-[960px] w-full mx-auto px-6 py-6">
		{@render children()}
	</main>

	<footer class="px-6 py-3 text-[11px] text-white/20 flex items-center justify-between">
		<span>SantaLucia 2004 · Mathews/Turner NN parameters</span>
		<span>Built by Vicenzo DeVito</span>
	</footer>
</div>
