type SharedValue = string | Blob | null | undefined;
type ShareKey = 'pdbcode' | 'biounit';

interface SharedCell {
  value: SharedValue;
  readonly: boolean;
  onchange: (value: SharedValue) => void;
  set: (value: string) => boolean;
}

interface ShareContext {
  self: {
    sharelist: string[];
    onsharelist: (keys: string[]) => void;
  };
  pdbcode?: SharedCell;
  biounit?: SharedCell;
}

type ShareAddress = string | number | null;
type ShareConnector = (
  updateServer?: ShareAddress,
  restServer?: ShareAddress,
  namespace?: string,
) => ShareContext;

declare global {
  interface Window {
    connect_seamless?: ShareConnector;
  }
}

export const seamlessState = $state({
  status: 'connecting' as 'connecting' | 'connected' | 'unavailable' | 'waiting-for-shares',
  ready: false,
  pdbCode: '',
  biounitIndex: 1,
  pdbCodePendingWrite: false,
  biounitIndexPendingWrite: false,
});

let client: ShareContext | undefined;
let pendingPdbCode: string | undefined;
let pendingBiounitIndex: number | undefined;
const received = new Set<ShareKey>();

function normalizePdbCode(value: SharedValue): string {
  return typeof value === 'string' ? value.trim().toUpperCase() : '';
}

function parseBiounitIndex(value: SharedValue): number {
  const raw = typeof value === 'string' ? value.trim() : '';
  const index = raw === '' ? Number.NaN : Number(raw);
  return Number.isSafeInteger(index) && index >= 1 ? index : 1;
}

function updateReady(): void {
  if (!received.has('pdbcode') || !received.has('biounit')) return;
  seamlessState.ready = true;
  seamlessState.status = 'connected';
}

function receivePdbCode(value: SharedValue): void {
  const code = normalizePdbCode(value);
  received.add('pdbcode');
  if (pendingPdbCode !== undefined && pendingPdbCode !== code) {
    seamlessState.pdbCode = pendingPdbCode;
  } else {
    pendingPdbCode = undefined;
    seamlessState.pdbCodePendingWrite = false;
    seamlessState.pdbCode = code;
  }
  updateReady();
}

function receiveBiounitIndex(value: SharedValue): void {
  const index = parseBiounitIndex(value);
  received.add('biounit');
  if (pendingBiounitIndex !== undefined && pendingBiounitIndex !== index) {
    seamlessState.biounitIndex = pendingBiounitIndex;
  } else {
    pendingBiounitIndex = undefined;
    seamlessState.biounitIndexPendingWrite = false;
    seamlessState.biounitIndex = index;
  }
  updateReady();
}

function attachShare(key: ShareKey, keys: string[]): void {
  if (!keys.includes(key)) return;
  const cell = client?.[key];
  if (!cell) return;
  cell.onchange = key === 'pdbcode' ? receivePdbCode : receiveBiounitIndex;
  if (cell.value !== undefined) cell.onchange(cell.value);
}

if (typeof window === 'undefined' || !window.connect_seamless) {
  seamlessState.status = 'unavailable';
} else {
  try {
    client = window.connect_seamless(null, null, 'preparic');
    client.self.onsharelist = keys => {
      attachShare('pdbcode', keys);
      attachShare('biounit', keys);
      if (!keys.includes('pdbcode') || !keys.includes('biounit')) {
        seamlessState.status = 'waiting-for-shares';
      }
    };
  } catch {
    seamlessState.status = 'unavailable';
  }
}

export function setSharedPdbCode(value: string): boolean {
  const code = value.trim().toUpperCase();
  const cell = client?.pdbcode;
  if (!client || !seamlessState.ready || !cell) {
    pendingPdbCode = code;
    seamlessState.pdbCodePendingWrite = true;
    seamlessState.pdbCode = code;
    return true;
  }
  if (cell.readonly) return false;
  if (normalizePdbCode(cell.value) === code) {
    pendingPdbCode = undefined;
    seamlessState.pdbCodePendingWrite = false;
    seamlessState.pdbCode = code;
    return true;
  }
  pendingPdbCode = code;
  seamlessState.pdbCodePendingWrite = true;
  const accepted = cell.set(code);
  if (accepted) seamlessState.pdbCode = code;
  return accepted;
}

export function setSharedBiounitIndex(value: number): boolean {
  const index = Number.isSafeInteger(value) && value >= 1 ? value : 1;
  const cell = client?.biounit;
  if (!client || !seamlessState.ready || !cell) {
    pendingBiounitIndex = index;
    seamlessState.biounitIndexPendingWrite = true;
    seamlessState.biounitIndex = index;
    return true;
  }
  if (cell.readonly) return false;
  if (parseBiounitIndex(cell.value) === index) {
    pendingBiounitIndex = undefined;
    seamlessState.biounitIndexPendingWrite = false;
    seamlessState.biounitIndex = index;
    return true;
  }
  pendingBiounitIndex = index;
  seamlessState.biounitIndexPendingWrite = true;
  const accepted = cell.set(String(index));
  if (accepted) seamlessState.biounitIndex = index;
  return accepted;
}
