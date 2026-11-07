export interface ExportBlock {
  code?: string;
  address?: string;
  message?: string;
  remedy?: string;
}

export interface ExportWarning {
  kind?: string;
  address?: string;
  message?: string;
}

export interface ExportAdmission {
  schema: 1;
  ready: boolean;
  blocks: ExportBlock[];
  warnings: ExportWarning[];
  digest: string | null;
}

export interface ExportReadinessInput {
  stage: number;
  connected: boolean;
  writePending: boolean;
  stage6EchoesFresh: boolean;
  admission: unknown;
  admissionState?: string;
  treeState?: string;
  admissionNodeChecksum?: string | null;
  admissionValueChecksum?: string | null;
  admissionError?: string;
  treeError?: string;
}

export interface ExportReadiness {
  enabled: boolean;
  reason: string;
  admission: ExportAdmission | null;
  warnings: ExportWarning[];
  requiresAcknowledgement: boolean;
}

function isRecord(value: unknown): value is Record<string, unknown> {
  return value !== null && typeof value === 'object' && !Array.isArray(value);
}

function readAdmission(value: unknown): ExportAdmission | null {
  if (!isRecord(value) || value.schema !== 1 || typeof value.ready !== 'boolean'
    || !Array.isArray(value.blocks) || !Array.isArray(value.warnings)
    || !(typeof value.digest === 'string' || value.digest === null)) return null;
  return value as unknown as ExportAdmission;
}

function blockReason(blocks: ExportBlock[]): string {
  return blocks.map(block => [block.address, block.message, block.remedy]
    .filter((part): part is string => typeof part === 'string' && part.length > 0)
    .join(': ')).filter(Boolean).join(' ');
}

/** Apply the browser-side Stage 6 gate without recalculating backend warnings. */
export function exportReadiness(input: ExportReadinessInput): ExportReadiness {
  const admission = readAdmission(input.admission);
  const warnings = admission?.warnings ?? [];
  const result = (reason: string, enabled = false): ExportReadiness => ({
    enabled,
    reason,
    admission,
    warnings,
    requiresAcknowledgement: warnings.length > 0,
  });

  if (input.stage !== 6) return result('Return to stage 6 to export.');
  if (!input.connected) return result('Reconnect to Preparic before exporting.');
  if (input.writePending) return result('Wait for pending project changes to be accepted.');
  if (!input.stage6EchoesFresh) return result('Wait for current stage 6 results.');

  if (input.admissionState !== 'complete') {
    if (input.admissionState === 'failed') {
      return result(`Export admission failed: ${input.admissionError || 'see the upstream error above.'}`);
    }
    if (input.admissionState === 'blocked') {
      return result(`An upstream calculation is blocking export${input.admissionError ? `: ${input.admissionError}` : '.'}`);
    }
    return result('Export admission is being recomputed.');
  }
  if (input.admissionNodeChecksum === null || input.admissionNodeChecksum === undefined
    || input.admissionValueChecksum === null || input.admissionValueChecksum === undefined
    || input.admissionNodeChecksum !== input.admissionValueChecksum) {
    return result('Waiting for the current export admission result.');
  }
  if (!admission) return result('The export admission result is malformed.');
  if (!admission.ready || admission.blocks.length > 0) {
    const detail = blockReason(admission.blocks);
    return result(detail ? `Export is blocked: ${detail}` : 'Export is blocked by admission checks.');
  }
  if (typeof admission.digest !== 'string' || admission.digest.length === 0) {
    return result('The export admission result has no snapshot digest.');
  }

  if (input.treeState === 'waiting' || input.treeState === 'computing' || !input.treeState) {
    return result('Export is being recomputed.');
  }
  if (input.treeState === 'blocked') {
    return result(`An upstream calculation is blocking export${input.treeError ? `: ${input.treeError}` : '.'}`);
  }
  if (input.treeState === 'failed') {
    return result(`Export tree failed: ${input.treeError || 'see the export error above.'}`);
  }
  if (input.treeState !== 'complete') return result('Waiting for the export tree to complete.');
  return result('', true);
}

export function exportNeedsAcknowledgement(admission: unknown): boolean {
  return (readAdmission(admission)?.warnings.length ?? 0) > 0;
}

/** The checkbox is checked only for the exact admission digest acknowledged in this tab. */
export function exportWarningAcknowledged(admission: unknown, acknowledgedDigest: string | null): boolean {
  const report = readAdmission(admission);
  if (!report) return false;
  if (report.warnings.length === 0) return true;
  return typeof report.digest === 'string' && report.digest.length > 0
    && acknowledgedDigest === report.digest;
}

export function exportIsStale(downloadedDigest: string | null, currentDigest: string | null): boolean {
  return typeof downloadedDigest === 'string' && downloadedDigest.length > 0
    && typeof currentDigest === 'string' && currentDigest.length > 0
    && downloadedDigest !== currentDigest;
}

function stableValue(value: unknown): unknown {
  if (Array.isArray(value)) return value.map(stableValue);
  if (!isRecord(value)) return value;
  return Object.fromEntries(Object.keys(value).sort().map(key => [key, stableValue(value[key])]));
}

/** A local cancellation key for decisions; it is not persisted or sent to Preparic. */
export function exportDecisionSignature(snapshot: unknown, energyThresholdEdits: unknown): string {
  return JSON.stringify(stableValue({ snapshot, energyThresholdEdits }));
}
