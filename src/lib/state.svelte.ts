import { cloneRangeRows, createRangeRow } from './ranges';
import type { ActiveCell, AssemblyInfo, ChainInfo, OtherChainsMode, RangeRow, Stage, Stacking } from './types';

export interface AppState {
  stage: Stage;
  pdbId: string;
  biounitIndex: number;
  assemblies: AssemblyInfo[];
  assemblyId: string;
  structureReady: boolean;
  loading: boolean;
  error: string;
  chains: ChainInfo[];
  rows: RangeRow[];
  activeCell: ActiveCell | null;
  otherChainsMode: OtherChainsMode;
  showSidechainSticks: boolean;
  stackings: Stacking[];
  stackingSelection: boolean[];
  staleStackingRevision: number | null;
  stage4OverlayFragments: number[];
  /** Stacking rows whose residues are hovered in the Mol* viewer. */
  viewerHoveredStackings: number[];
  /** Stacking row hovered in the table. */
  tableHoveredStacking: number | null;
}

export const appState = $state<AppState>({
  stage: 1,
  pdbId: '',
  biounitIndex: 1,
  assemblies: [],
  assemblyId: '',
  structureReady: false,
  loading: false,
  error: '',
  chains: [],
  rows: [],
  activeCell: null,
  otherChainsMode: 'transparent',
  showSidechainSticks: false,
  stackings: [],
  stackingSelection: [],
  staleStackingRevision: null,
  stage4OverlayFragments: [],
  viewerHoveredStackings: [],
  tableHoveredStacking: null,
});

export function enterStage2(initialRows?: RangeRow[]): void {
  appState.rows = initialRows?.length ? cloneRangeRows(initialRows) : [createRangeRow()];
  appState.activeCell = null;
  appState.showSidechainSticks = false;
  appState.stage = 2;
}

export function enterStage3(): void {
  appState.activeCell = null;
  appState.stage = 3;
}

export function returnToStage2(): void {
  appState.activeCell = null;
  appState.stage = 2;
}

export function clearStackingState(): void {
  appState.stackings = [];
  appState.stackingSelection = [];
  appState.staleStackingRevision = null;
}

export function invalidateStackingState(staleRevision: number): void {
  clearStackingState();
  appState.staleStackingRevision = staleRevision;
}

export function clearStructureDependentState(): void {
  appState.rows = [];
  appState.activeCell = null;
  appState.chains = [];
  appState.showSidechainSticks = false;
  clearStackingState();
}

export function returnToStage1(): void {
  clearStructureDependentState();
  appState.stage = 1;
}
