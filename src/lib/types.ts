export type Stage = 1 | 2 | 3;

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

export type Vec3 = [number, number, number];

export interface Stacking {
  rna_resname: string;
  rna_resid: number;
  protein_range_pdb_id: number;
  protein_resid: number;
  distance: number;
  tilt_angle: number;
  slip_angle: number;
  rna_centroid: Vec3;
  protein_centroid: Vec3;
}

/** Resolved participant geometry passed to the viewer by the stage controller. */
export interface StackingViewerInput {
  rna: RangeSpec;
  protein: RangeSpec;
  rnaCentroid: Vec3;
  proteinCentroid: Vec3;
}

export interface RenderSpec {
  focus: ChainKey | null;
  range?: RangeSpec;
  alwaysVisible?: RangeSpec[];
  sidechainRanges?: RangeSpec[];
  stackings?: StackingViewerInput[];
  others?: OtherChainsMode;
}
