<script>
	import '../app.css';
	import { onMount } from 'svelte';

	let { children } = $props();
	let theme = $state('cream');

	onMount(() => {
		const saved = localStorage.getItem('pw-theme');
		theme = saved === 'dark' ? 'dark' : 'cream';
		document.documentElement.dataset.theme = theme;
	});

	function toggleTheme() {
		theme = theme === 'dark' ? 'cream' : 'dark';
		document.documentElement.dataset.theme = theme;
		localStorage.setItem('pw-theme', theme);
	}
</script>

<div class="min-h-screen flex flex-col">
	<div class="mx-auto w-full" style="max-width: 1180px; padding: var(--pad-page-y) var(--pad-page-x);">
		<header class="flex items-center justify-between" style="margin-bottom: 28px;">
			<div class="pw-eyebrow" style="font-size: 12px; letter-spacing: 0.12em;">
				Primer Workbench
			</div>
			<div class="flex items-center gap-3">
				<button class="pw-btn-quiet" onclick={toggleTheme} title="Toggle theme">
					{theme === 'dark' ? 'Light' : 'Dark'}
				</button>
				<a href="https://github.com/V-DeVito/PrimerDesigner" target="_blank" rel="noopener"
				   class="pw-btn-quiet" style="text-decoration: none;">
					GitHub
				</a>
			</div>
		</header>

		<main>
			{@render children()}
		</main>
	</div>
</div>
