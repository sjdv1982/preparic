<script lang="ts">
  import { tick } from 'svelte';
  import { appState } from '../lib/state.svelte';
  import { addRow, isCellInvalid, moveCell, moveRow, removeRow, setActiveCell, updateRow, validateRows } from '../lib/ranges';
  import type { ActiveCell, RangeColumn } from '../lib/types';
  import type { ValidationIssue } from '../lib/ranges';

  let { issues }: { issues: ValidationIssue[] } = $props();
  let editing: ActiveCell | null = $state(null);

  function activate(next: ActiveCell | null) {
    appState.rows = setActiveCell(appState.rows, appState.activeCell, next);
    appState.activeCell = next;
    if (editing && (!next || editing.rowId !== next.rowId || editing.column !== next.column)) editing = null;
  }

  function cellValue(rowId: string, column: RangeColumn): string {
    const row = appState.rows.find(item => item.id === rowId);
    return row?.[column] ?? '';
  }

  async function focusCell(cell: ActiveCell) {
    activate(cell);
    await tick();
    const element = document.querySelector<HTMLElement>(`[data-cell-key="${cell.rowId}:${cell.column}"]`);
    element?.focus();
  }

  async function beginEdit(cell: ActiveCell, firstKey?: string) {
    activate(cell);
    editing = cell;
    if (firstKey && cell.column !== 'chain') {
      appState.rows = updateRow(appState.rows, cell.rowId, { [cell.column]: firstKey });
    }
    await tick();
    const element = document.querySelector<HTMLInputElement | HTMLSelectElement>(`[data-editor-key="${cell.rowId}:${cell.column}"]`);
    element?.focus();
    if (element instanceof HTMLInputElement && !firstKey) element.select();
  }

  function updateValue(rowId: string, column: RangeColumn, value: string) {
    appState.rows = updateRow(appState.rows, rowId, { [column]: value });
  }

  function clearCell(rowId: string, column: RangeColumn) {
    updateValue(rowId, column, '');
    editing = null;
  }

  function moveActive(cell: ActiveCell, key: 'ArrowUp' | 'ArrowDown' | 'ArrowLeft' | 'ArrowRight' | 'Tab', shift = false) {
    const next = moveCell(appState.rows, cell, key, shift);
    if (next) void focusCell(next);
  }

  function handleCellKeydown(event: KeyboardEvent, cell: ActiveCell) {
    if (event.key === 'Enter') {
      event.preventDefault();
      editing = null;
      return;
    }
    if (event.key === 'Tab' || event.key.startsWith('Arrow')) {
      event.preventDefault();
      editing = null;
      moveActive(cell, event.key as 'ArrowUp' | 'ArrowDown' | 'ArrowLeft' | 'ArrowRight' | 'Tab', event.shiftKey);
      return;
    }
    if (event.key === 'Delete' || event.key === 'Backspace') {
      event.preventDefault();
      clearCell(cell.rowId, cell.column);
      return;
    }
    if (event.key.length === 1 && !event.ctrlKey && !event.metaKey && !event.altKey) {
      event.preventDefault();
      if (cell.column === 'chain') {
        const match = appState.chains.find(chain => chain.name.toLowerCase().startsWith(event.key.toLowerCase()));
        if (match) updateValue(cell.rowId, 'chain', match.name);
        return;
      }
      if (/^[0-9-]$/.test(event.key)) void beginEdit(cell, event.key);
    }
  }

  function handleEditorKeydown(event: KeyboardEvent, cell: ActiveCell) {
    if (event.key === 'Escape') {
      event.preventDefault();
      editing = null;
      void focusCell(cell);
    } else {
      handleCellKeydown(event, cell);
    }
  }

  function toggleAlwaysVisible(rowId: string, checked: boolean) {
    appState.rows = updateRow(appState.rows, rowId, { alwaysVisible: checked });
  }

  function clearRange(rowId: string) {
    appState.rows = updateRow(appState.rows, rowId, { start: '', end: '', alwaysVisible: false });
    editing = null;
  }

  function addBelowActive() {
    const activeRowId = appState.activeCell?.rowId;
    appState.rows = addRow(appState.rows, activeRowId);
  }

  function deleteRow(rowId: string) {
    appState.rows = removeRow(appState.rows, rowId);
    if (appState.activeCell?.rowId === rowId) activate(null);
  }

  function handleOutsidePointerDown(event: PointerEvent) {
    const target = event.target;
    if (!(target instanceof Element)) return;
    if (target.closest('.range-table') || target.closest('.viewer-shell') || target.closest('.table-actions') || target.closest('.sidechain-toggle')) return;
    activate(null);
  }

  const chainNames = $derived(appState.chains.map(chain => chain.name));
  const rowIssues = $derived(validateRows(appState.rows, appState.chains, true));
</script>

<svelte:window onpointerdown={handleOutsidePointerDown} />

<div class="range-table-wrap">
  <table class="range-table" role="grid" aria-label="Residue ranges">
    <thead>
      <tr>
        <th class="pin-heading" title="Always visible">Pin</th>
        <th>Chain</th>
        <th>Start</th>
        <th>End</th>
        <th class="controls-heading"><span class="sr-only">Row controls</span></th>
      </tr>
    </thead>
    <tbody>
      {#each appState.rows as row, rowIndex (row.id)}
        <tr>
          <td class="pin-cell">
            <input type="checkbox" aria-label={`Keep range ${rowIndex + 1} always visible`} checked={row.alwaysVisible} onchange={(event) => toggleAlwaysVisible(row.id, event.currentTarget.checked)} />
          </td>
          {#each ['chain', 'start', 'end'] as columnValue (columnValue)}
            {@const column = columnValue as RangeColumn}
            {@const cell = { rowId: row.id, column }}
            {@const active = appState.activeCell?.rowId === row.id && appState.activeCell.column === column}
            {@const editingThis = editing?.rowId === row.id && editing.column === column}
            {@const invalid = isCellInvalid(issues, row.id, column)}
            <td class:active-cell={active} class:invalid-cell={invalid} title={rowIssues.find(issue => issue.rowId === row.id && issue.column === column)?.message ?? ''}>
              {#if editingThis && column === 'chain'}
                <select class="cell-editor" data-editor-key={`${row.id}:${column}`} aria-label={`Chain for row ${rowIndex + 1}`} value={row.chain} onchange={(event) => { updateValue(row.id, column, event.currentTarget.value); editing = null; }} onkeydown={(event) => handleEditorKeydown(event, cell)} onblur={() => { editing = null; }}>
                  <option value="">Choose…</option>
                  {#each chainNames as name (name)}<option value={name}>{name}</option>{/each}
                </select>
              {:else if editingThis}
                <input class="cell-editor" data-editor-key={`${row.id}:${column}`} aria-label={`${column === 'start' ? 'Start' : 'End'} residue for row ${rowIndex + 1}`} inputmode="numeric" value={row[column]} oninput={(event) => updateValue(row.id, column, event.currentTarget.value)} onkeydown={(event) => handleEditorKeydown(event, cell)} onblur={() => { editing = null; }} />
              {:else}
                <button type="button" class="cell-display" class:cell-active={active} data-cell-key={`${row.id}:${column}`} tabindex={active || !appState.activeCell ? 0 : -1} aria-label={`${column} ${row[column] || 'empty'}, row ${rowIndex + 1}`} onclick={() => activate(cell)} ondblclick={() => void beginEdit(cell)} onkeydown={(event) => handleCellKeydown(event, cell)}>
                  {#if row[column]}
                    {row[column]}
                  {:else}
                    <span class="cell-placeholder">{column === 'chain' ? 'Choose…' : '—'}</span>
                  {/if}
                </button>
              {/if}
            </td>
          {/each}
          <td class="row-controls">
            <button type="button" class="icon-button" aria-label={`Move row ${rowIndex + 1} up`} title="Move up" disabled={rowIndex === 0} onclick={() => appState.rows = moveRow(appState.rows, row.id, -1)}>↑</button>
            <button type="button" class="icon-button" aria-label={`Move row ${rowIndex + 1} down`} title="Move down" disabled={rowIndex === appState.rows.length - 1} onclick={() => appState.rows = moveRow(appState.rows, row.id, 1)}>↓</button>
            <button type="button" class="icon-button clear-range-button" aria-label={`Clear range for row ${rowIndex + 1}`} title="Clear range" onclick={() => clearRange(row.id)}>↺</button>
            <button type="button" class="icon-button remove-button" aria-label={`Remove row ${rowIndex + 1}`} title="Remove row" onclick={() => deleteRow(row.id)}>×</button>
          </td>
        </tr>
      {/each}
    </tbody>
  </table>
  <div class="table-actions">
    <button type="button" class="add-row-button" onclick={addBelowActive}><span aria-hidden="true">+</span> Add row</button>
    <span class="table-key-help">Click to select <span>·</span> Double-click or type to edit</span>
  </div>
</div>
