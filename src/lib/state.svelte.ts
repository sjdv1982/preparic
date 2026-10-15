import { createRangeRow } from './ranges';
import type { ActiveCell, AssemblyInfo, ChainInfo, OtherChainsMode, RangeRow, Stage } from './types';

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
});

export function enterStage2(): void {
  appState.rows = [createRangeRow()];
  appState.activeCell = null;
  appState.showSidechainSticks = false;
  appState.stage = 2;
}

export function returnToStage1(): void {
  appState.rows = [];
  appState.activeCell = null;
  appState.chains = [];
  appState.showSidechainSticks = false;
  appState.stage = 1;
}
