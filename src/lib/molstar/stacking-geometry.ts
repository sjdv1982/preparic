import type { Vec3 } from '../types';

export interface DashSegment {
  start: Vec3;
  end: Vec3;
}

/** Split a centroid connector into short cylinders separated by visible gaps. */
export function stackingDashSegments(
  start: Vec3,
  end: Vec3,
  dashCount = 8,
  dashFraction = 0.58,
): DashSegment[] {
  if (!Number.isSafeInteger(dashCount) || dashCount <= 0) return [];
  if (!Number.isFinite(dashFraction) || dashFraction <= 0 || dashFraction > 1) return [];
  if (![...start, ...end].every(Number.isFinite)) return [];

  const delta: Vec3 = [end[0] - start[0], end[1] - start[1], end[2] - start[2]];
  const lengthSquared = delta[0] ** 2 + delta[1] ** 2 + delta[2] ** 2;
  if (!Number.isFinite(lengthSquared) || lengthSquared <= 0) return [];

  return Array.from({ length: dashCount }, (_, index) => {
    const from = index / dashCount;
    const to = from + dashFraction / dashCount;
    return {
      start: [
        start[0] + delta[0] * from,
        start[1] + delta[1] * from,
        start[2] + delta[2] * from,
      ],
      end: [
        start[0] + delta[0] * to,
        start[1] + delta[1] * to,
        start[2] + delta[2] * to,
      ],
    };
  });
}
