<script lang="ts">
  import { onMount } from 'svelte';
  import type { MolViewer } from '../lib/molstar/viewer';
  import type { HoverInfo, PickInfo } from '../lib/types';

  let { onReady, onError, onPick, onHover, chainName }: {
    onReady: (viewer: MolViewer) => void;
    onError: (message: string) => void;
    onPick: (pick: PickInfo) => void;
    onHover?: (hover: HoverInfo | null) => void;
    chainName?: (hover: HoverInfo) => string;
  } = $props();

  let container: HTMLDivElement;
  let canvas: HTMLCanvasElement;
  let viewer = $state<MolViewer | null>(null);

  let hover = $state<HoverInfo | null>(null);
  let pointer = $state({ x: 0, y: 0 });

  $effect(() => {
    viewer?.setOnPick(onPick);
  });

  $effect(() => {
    viewer?.setOnHover(value => {
      hover = value;
      onHover?.(value);
    });
  });

  function trackPointer(event: PointerEvent) {
    const rect = container.getBoundingClientRect();
    pointer = { x: event.clientX - rect.left, y: event.clientY - rect.top };
  }

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

<div class="viewer-shell" role="presentation" bind:this={container} onpointermove={trackPointer} onpointerleave={() => { hover = null; onHover?.(null); }}>
  <canvas bind:this={canvas} aria-label="Interactive molecular structure viewer"></canvas>
  <div class="viewer-hint"><span class="hint-dot"></span>Drag to rotate <span class="hint-separator">·</span> Scroll to zoom <span class="hint-separator">·</span> Click to pick</div>
  {#if hover}
    <div class="viewer-tooltip" role="tooltip" style:left="{pointer.x}px" style:top="{pointer.y}px">{chainName?.(hover) ?? hover.labelAsymId} · {hover.resname} {hover.labelSeqId}</div>
  {/if}
</div>
