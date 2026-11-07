import { describe, expect, it, vi } from 'vitest';
import {
  cancelExportRequest,
  createExportRequest,
  exportRequestIsCurrent,
  saveExportBlob,
  transferExport,
} from '../../src/lib/export-download';

const digest = 'a'.repeat(64);
const zipBytes = Uint8Array.from([0x50, 0x4b, 0x03, 0x04]);

function response(status = 200): Response {
  return new Response(status === 200 ? zipBytes : 'not ready', { status });
}

function transferOptions(overrides: Record<string, any> = {}) {
  const request = createExportRequest(digest);
  const fetcher = vi.fn(async (_url: any, _init: any) => response());
  const readManifest = vi.fn(async () => ({ schema: 1, snapshot_sha256: digest }));
  const save = vi.fn();
  return {
    request,
    fetcher,
    readManifest,
    save,
    options: {
      url: '/preparic/export.zip',
      digest,
      filename: '1b7f-alaric-aaaaaaaa.zip',
      request,
      fetcher,
      readManifest,
      isCurrent: () => exportRequestIsCurrent(request, request, digest, true),
      save,
      ...overrides,
    },
  };
}

describe('E6 export request state', () => {
  it('allows only the active uncancelled request for the current ready digest', () => {
    const request = createExportRequest(digest);
    expect(exportRequestIsCurrent(request, request, digest, true)).toBe(true);
    expect(exportRequestIsCurrent(request, null, digest, true)).toBe(false);
    expect(exportRequestIsCurrent(request, request, 'd'.repeat(64), true)).toBe(false);
    expect(exportRequestIsCurrent(request, request, digest, false)).toBe(false);
  });

  it('aborts and permanently cancels a request, even if the digest later reverts', () => {
    const request = createExportRequest(digest);
    expect(exportRequestIsCurrent(request, request, 'd'.repeat(64), true)).toBe(false);
    cancelExportRequest(request);
    expect(request.cancelled).toBe(true);
    expect(request.controller.signal.aborted).toBe(true);
    expect(exportRequestIsCurrent(request, request, digest, true)).toBe(false);
  });
});

describe('E6 guarded export transfer', () => {
  it('fetches without consent, disables cache, validates the manifest, then saves the named archive', async () => {
    const { options, fetcher, readManifest, save } = transferOptions();
    const result = await transferExport(options);

    expect(result.kind).toBe('downloaded');
    expect(fetcher).toHaveBeenCalledOnce();
    const init = fetcher.mock.calls[0][1];
    expect(init.cache).toBe('no-store');
    expect(init.signal).toBe(options.request.controller.signal);
    expect(Array.from(new Headers(init.headers).keys())).toEqual([]);
    expect(readManifest).toHaveBeenCalledOnce();
    expect(save).toHaveBeenCalledOnce();
    expect(save).toHaveBeenCalledWith(expect.anything(), '1b7f-alaric-aaaaaaaa.zip');
  });

  it.each([404, 409, 503])('surfaces actionable HTTP %s without releasing an archive', async status => {
    const fetcher = vi.fn(async (_url: any, _init: any) => response(status));
    const { options, save, readManifest } = transferOptions({ fetcher });
    const result = await transferExport(options);

    expect(result.kind).toBe('error');
    expect(result.message).toContain(String(status));
    expect(result.message.length).toBeGreaterThan(0);
    expect(readManifest).not.toHaveBeenCalled();
    expect(save).not.toHaveBeenCalled();
  });

  it('rejects other non-200 statuses without attempting ZIP parsing or saving', async () => {
    const fetcher = vi.fn(async (_url: any, _init: any) => response(418));
    const { options, save, readManifest } = transferOptions({ fetcher });
    const result = await transferExport(options);
    expect(result.kind).toBe('error');
    expect(result.message).toContain('418');
    expect(readManifest).not.toHaveBeenCalled();
    expect(save).not.toHaveBeenCalled();
  });

  it('rejects an invalid ZIP or malformed manifest without saving', async () => {
    const { options, save } = transferOptions({
      readManifest: vi.fn(async () => { throw new Error('Invalid ZIP archive.'); }),
    });
    const result = await transferExport(options);
    expect(result.kind).toBe('error');
    expect(save).not.toHaveBeenCalled();

    const malformed = transferOptions({ readManifest: vi.fn(async () => ({ schema: 2, snapshot_sha256: digest })) });
    const malformedResult = await transferExport(malformed.options);
    expect(malformedResult.kind).toBe('error');
    expect(malformed.save).not.toHaveBeenCalled();
  });

  it('rejects a manifest for a different digest and a digest that became stale during transfer', async () => {
    const mismatch = transferOptions({
      readManifest: vi.fn(async () => ({ schema: 1, snapshot_sha256: 'd'.repeat(64) })),
    });
    expect((await transferExport(mismatch.options)).kind).toBe('stale');
    expect(mismatch.save).not.toHaveBeenCalled();

    const changed = transferOptions({ isCurrent: () => false });
    const changedResult = await transferExport(changed.options);
    expect(changedResult.kind).toBe('stale');
    expect(changedResult.message).toContain('Project changed during export');
    expect(changed.save).not.toHaveBeenCalled();
  });

  it('discards a late 200 response after cancellation and forwards the abort signal', async () => {
    let resolveResponse!: (value: Response) => void;
    const pendingResponse = new Promise<Response>(resolve => { resolveResponse = resolve; });
    const fetcher = vi.fn((_url: any, _init: any) => pendingResponse);
    const { options, request, save } = transferOptions({ fetcher });
    const pending = transferExport(options);
    expect(fetcher.mock.calls[0][1].signal).toBe(request.controller.signal);

    cancelExportRequest(request);
    resolveResponse(response());
    expect((await pending).kind).toBe('cancelled');
    expect(save).not.toHaveBeenCalled();
  });
});

describe('E6 ZIP save cleanup', () => {
  it('removes its temporary anchor and revokes its object URL after the click', () => {
    const anchor = { href: '', download: '', click: vi.fn() };
    const createObjectURL = vi.fn(() => 'blob:export');
    const revokeObjectURL = vi.fn();
    const appendAnchor = vi.fn();
    const removeAnchor = vi.fn();

    saveExportBlob(new Blob([zipBytes], { type: 'application/zip' }), 'export.zip', {
      createObjectURL,
      revokeObjectURL,
      createAnchor: () => anchor,
      appendAnchor,
      removeAnchor,
    });

    expect(anchor.href).toBe('blob:export');
    expect(anchor.download).toBe('export.zip');
    expect(appendAnchor).toHaveBeenCalledWith(anchor);
    expect(anchor.click).toHaveBeenCalledOnce();
    expect(removeAnchor).toHaveBeenCalledWith(anchor);
    expect(revokeObjectURL).toHaveBeenCalledWith('blob:export');
  });

  it('still removes the anchor and revokes the URL when the browser click fails', () => {
    const anchor = { href: '', download: '', click: vi.fn(() => { throw new Error('click failed'); }) };
    const revokeObjectURL = vi.fn();
    const removeAnchor = vi.fn();

    expect(() => saveExportBlob(new Blob([zipBytes]), 'export.zip', {
      createObjectURL: () => 'blob:export',
      revokeObjectURL,
      createAnchor: () => anchor,
      appendAnchor: vi.fn(),
      removeAnchor,
    })).toThrow('click failed');

    expect(removeAnchor).toHaveBeenCalledWith(anchor);
    expect(revokeObjectURL).toHaveBeenCalledWith('blob:export');
  });

  it('revokes the URL even if removing the temporary anchor throws', () => {
    const anchor = { href: '', download: '', click: vi.fn() };
    const revokeObjectURL = vi.fn();

    expect(() => saveExportBlob(new Blob([zipBytes]), 'export.zip', {
      createObjectURL: () => 'blob:export',
      revokeObjectURL,
      createAnchor: () => anchor,
      appendAnchor: vi.fn(),
      removeAnchor: vi.fn(() => { throw new Error('remove failed'); }),
    })).toThrow('remove failed');

    expect(revokeObjectURL).toHaveBeenCalledWith('blob:export');
  });
});
