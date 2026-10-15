export type Stage = 1 | 2;

export interface ChainKey {
  labelAsymId: string;
  instanceId: string;
}

export interface ChainInfo extends ChainKey {
  name: string;
  seqIds: Set<number>;
  seqMin: number;
  seqMax: number;
}

export type RangeColumn = 'chain' | 'start' | 'end';

export interface ActiveCell {
  rowId: string;
  column: RangeColumn;
}

export interface RangeRow {
  id: string;
  alwaysVisible: boolean;
  chain: string;
  start: string;
  end: string;
}

export type OtherChainsMode = 'visible' | 'transparent' | 'invisible';

export interface AssemblyParticipant {
  labelAsymId: string;
  count: number;
}

export interface AssemblyInfo {
  id: string;
  label: string;
  details?: string;
  participants: AssemblyParticipant[];
}

export interface PickInfo extends ChainKey {
  labelSeqId: number;
}

export interface RangeSpec extends ChainKey {
  start: number;
  end: number;
}

export interface RenderSpec {
  focus: ChainKey | null;
  range?: RangeSpec;
  alwaysVisible?: RangeSpec[];
  sidechainRanges?: RangeSpec[];
  others?: OtherChainsMode;
}
