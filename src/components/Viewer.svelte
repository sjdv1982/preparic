<script lang="ts">
  import { onMount } from 'svelte';
  import type { MolViewer } from '../lib/molstar/viewer';
  import type { PickInfo } from '../lib/types';

  let { onReady, onError, onPick }: {
    onReady: (viewer: MolViewer) => void;
    onError: (message: string) => void;
    onPick: (pick: PickInfo) => void;
  } = $props();

  let container: HTMLDivElement;
  let canvas: HTMLCanvasElement;
  let viewer = $state<MolViewer | null>(null);

  $effect(() => {
    viewer?.setOnPick(onPick);
  });

  onMount(() => {
    let disposed = false;
    void import('../lib/molstar/viewer').then(({ MolViewer }) => MolViewer.create(canvas, container)).then(created => {
      if (disposed) {
        created.dispose();
        return;
      }
      viewer = created;
      onReady(created);
    }).catch(error => {
      if (!disposed) onError(error instanceof Error ? error.message : String(error));
    });

    return () => {
      disposed = true;
      viewer?.dispose();
    };
  });
</script>

<div class="viewer-shell" bind:this={container}>
  <canvas bind:this={canvas} aria-label="Interactive molecular structure viewer"></canvas>
  <div class="viewer-hint"><span class="hint-dot"></span>Drag to rotate <span class="hint-separator">·</span> Scroll to zoom <span class="hint-separator">·</span> Click to pick</div>
</div>
