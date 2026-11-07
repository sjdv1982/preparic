import { unzipSync } from 'fflate';

export interface ExportRequest {
  digest: string;
  controller: AbortController;
  cancelled: boolean;
}

export type ExportTransferOutcome =
  | { kind: 'downloaded'; message: string }
  | { kind: 'cancelled'; message: string }
  | { kind: 'stale'; message: string }
  | { kind: 'error'; message: string; status?: number };

export interface ExportTransferOptions {
  url: string;
  digest: string;
  filename: string;
  request: ExportRequest;
  fetcher?: (input: string, init: RequestInit) => Promise<Response>;
  readManifest?: (zip: Blob) => Promise<unknown>;
  isCurrent: () => boolean;
  save: (zip: Blob, filename: string) => void;
}

export interface ExportDownloadEnvironment {
  createObjectURL: (blob: Blob) => string;
  revokeObjectURL: (url: string) => void;
  createAnchor: () => { href: string; download: string; click: () => void };
  appendAnchor: (anchor: { href: string; download: string; click: () => void }) => void;
  removeAnchor: (anchor: { href: string; download: string; click: () => void }) => void;
}

export function createExportRequest(digest: string, controller = new AbortController()): ExportRequest {
  return { digest, controller, cancelled: false };
}

export function cancelExportRequest(request: ExportRequest): void {
  if (request.cancelled) return;
  request.cancelled = true;
  request.controller.abort();
}

export function exportRequestIsCurrent(
  request: ExportRequest,
  activeRequest: ExportRequest | null,
  currentDigest: string | null,
  stillReady: boolean,
): boolean {
  return !request.cancelled && activeRequest === request && request.digest === currentDigest && stillReady;
}

async function defaultReadManifest(zip: Blob): Promise<unknown> {
  const bytes = new Uint8Array(await zip.arrayBuffer());
  const files = unzipSync(bytes);
  const manifest = files['alaric_project/DATA/export-manifest.json'];
  if (!manifest) throw new Error('The ZIP is missing alaric_project/DATA/export-manifest.json.');
  try {
    return JSON.parse(new TextDecoder().decode(manifest)) as unknown;
  } catch {
    throw new Error('The export manifest is not valid JSON.');
  }
}

function isSchemaOneManifest(value: unknown): value is { schema: 1; snapshot_sha256: unknown } {
  return value !== null && typeof value === 'object' && !Array.isArray(value)
    && (value as Record<string, unknown>).schema === 1
    && Object.hasOwn(value, 'snapshot_sha256');
}

async function readErrorBody(response: Response): Promise<unknown> {
  try { return await response.json(); } catch { return null; }
}

function treeFailureReason(value: unknown): { state: string; message: string } {
  if (!value || typeof value !== 'object' || Array.isArray(value)) {
    return { state: typeof value === 'string' ? value : '', message: '' };
  }
  const record = value as Record<string, unknown>;
  const state = typeof record.state === 'string' ? record.state : '';
  const message = typeof record.exception === 'string' ? record.exception
    : typeof record.message === 'string' ? record.message : '';
  return { state, message };
}

async function errorForResponse(response: Response): Promise<ExportTransferOutcome> {
  const status = response.status;
  if (status === 409) {
    const body = await readErrorBody(response);
    const cells = body && typeof body === 'object' && !Array.isArray(body)
      ? (body as Record<string, unknown>).cells : null;
    const tree = cells && typeof cells === 'object' && !Array.isArray(cells)
      ? (cells as Record<string, unknown>).export_tree : null;
    const { state, message } = treeFailureReason(tree);
    if (state === 'waiting' || state === 'computing') {
      return { kind: 'error', status, message: `Export is being recomputed (HTTP ${status}). Wait for it to finish, then try again.` };
    }
    if (state === 'blocked') {
      return { kind: 'error', status, message: `An upstream calculation is blocking export (HTTP ${status}). Resolve the displayed Preparic error, then try again.` };
    }
    if (state === 'failed') {
      return { kind: 'error', status, message: `The export tree failed (HTTP ${status})${message ? `: ${message}` : '. Inspect the export error above before trying again.'}` };
    }
    return { kind: 'error', status, message: `Export is not ready (HTTP ${status}). Review the current export state, then try again.` };
  }
  if (status === 503) {
    return { kind: 'error', status, message: 'Preparic could not resolve the export archive (HTTP 503). Check the connection and try again.' };
  }
  if (status === 404) {
    return { kind: 'error', status, message: 'The export download is not registered (HTTP 404). Wait for the Preparic backend to restart, then try again.' };
  }
  return { kind: 'error', status, message: `Export download returned HTTP ${status}. No archive was saved.` };
}

function staleOutcome(): ExportTransferOutcome {
  return { kind: 'stale', message: 'Project changed during export; export again.' };
}

/** Fetch, inspect the embedded manifest, and release a ZIP only for the same live request. */
export async function transferExport(options: ExportTransferOptions): Promise<ExportTransferOutcome> {
  if (options.request.cancelled) return { kind: 'cancelled', message: '' };
  if (!options.isCurrent()) return staleOutcome();

  const fetcher = options.fetcher ?? fetch;
  try {
    const response = await fetcher(options.url, {
      cache: 'no-store',
      signal: options.request.controller.signal,
    });
    if (options.request.cancelled) return { kind: 'cancelled', message: '' };
    if (!options.isCurrent()) return staleOutcome();
    if (response.status !== 200) return await errorForResponse(response);

    const responseBlob = await response.blob();
    let manifest: unknown;
    try {
      manifest = await (options.readManifest ?? defaultReadManifest)(responseBlob);
    } catch (error) {
      const detail = error instanceof Error ? error.message : 'The ZIP could not be read.';
      return { kind: 'error', status: 200, message: `The export ZIP is invalid: ${detail}` };
    }
    if (!isSchemaOneManifest(manifest)) {
      return { kind: 'error', status: 200, message: 'The export ZIP has a malformed schema-1 manifest. No archive was saved.' };
    }
    if (manifest.snapshot_sha256 !== options.digest) return staleOutcome();
    if (options.request.cancelled) return { kind: 'cancelled', message: '' };
    if (!options.isCurrent()) return staleOutcome();

    const zip = new Blob([responseBlob], { type: 'application/zip' });
    options.save(zip, options.filename);
    return { kind: 'downloaded', message: '' };
  } catch (error) {
    if (options.request.cancelled || (error instanceof Error && error.name === 'AbortError')) {
      return { kind: 'cancelled', message: '' };
    }
    const detail = error instanceof Error ? error.message : 'The transfer could not be completed.';
    return { kind: 'error', message: `Export failed: ${detail}` };
  }
}

function browserEnvironment(): ExportDownloadEnvironment {
  return {
    createObjectURL: blob => URL.createObjectURL(blob),
    revokeObjectURL: url => URL.revokeObjectURL(url),
    createAnchor: () => document.createElement('a'),
    appendAnchor: anchor => document.body.appendChild(anchor as HTMLAnchorElement),
    removeAnchor: anchor => {
      const element = anchor as HTMLAnchorElement;
      if (document.body.contains(element)) document.body.removeChild(element);
    },
  };
}

/** Trigger a browser download and release its temporary DOM and object URL even on failure. */
export function saveExportBlob(
  blob: Blob,
  filename: string,
  environment: ExportDownloadEnvironment = browserEnvironment(),
): void {
  let objectUrl: string | undefined;
  let anchor: ReturnType<ExportDownloadEnvironment['createAnchor']> | undefined;
  let appended = false;
  try {
    objectUrl = environment.createObjectURL(blob);
    anchor = environment.createAnchor();
    anchor.href = objectUrl;
    anchor.download = filename;
    environment.appendAnchor(anchor);
    appended = true;
    anchor.click();
  } finally {
    try {
      if (appended && anchor) environment.removeAnchor(anchor);
    } finally {
      if (objectUrl !== undefined) environment.revokeObjectURL(objectUrl);
    }
  }
}
