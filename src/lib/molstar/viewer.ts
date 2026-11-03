import { StructureElement, StructureProperties, Unit } from 'molstar/lib/mol-model/structure.js';
import { MmcifFormat } from 'molstar/lib/mol-model-formats/structure/mmcif.js';
import { ModelSymmetry } from 'molstar/lib/mol-model-formats/structure/property/symmetry.js';
import { PluginBehaviors } from 'molstar/lib/mol-plugin/behavior.js';
import { PluginContext } from 'molstar/lib/mol-plugin/context.js';
import { DefaultPluginSpec } from 'molstar/lib/mol-plugin/spec.js';
import { Mesh } from 'molstar/lib/mol-geo/geometry/mesh/mesh.js';
import { MeshBuilder } from 'molstar/lib/mol-geo/geometry/mesh/mesh-builder.js';
import { addSimpleCylinder } from 'molstar/lib/mol-geo/geometry/mesh/builder/cylinder.js';
import { Shape } from 'molstar/lib/mol-model/shape.js';
import { Vec3 as MolVec3 } from 'molstar/lib/mol-math/linear-algebra.js';
import { Color } from 'molstar/lib/mol-util/color/color.js';
import { ParamDefinition as PD } from 'molstar/lib/mol-util/param-definition.js';
import { StateTransformer } from 'molstar/lib/mol-state/index.js';
import { PluginStateObject as SO } from 'molstar/lib/mol-plugin-state/objects.js';
import { ShapeRepresentation3D } from 'molstar/lib/mol-plugin-state/transforms/representation.js';
import { StateSelection } from 'molstar/lib/mol-state/index.js';
import type { ResidueIndex } from 'molstar/lib/mol-model/structure/model/indexing.js';
import { chainExpression, complementExpression, rangeExpression, sidechainRangeExpression, unionExpression } from './expressions';
import { stackingDashSegments, type DashSegment } from './stacking-geometry';
import { nameChains, type ChainUnitInput } from '../chains';
import { uniqueStackingParticipantRanges } from '../stackings';
import type { AssemblyInfo, ChainInfo, PickInfo, RenderSpec, StackingViewerInput, Vec3 } from '../types';

const RENDER_TAG = 'preparic-phase1-render';
const PREPARED_PDB_TAG = 'preparic-prepared-pdb';
interface StructureSource {
  cif: string;
  id: string;
  dataRef: any;
  trajectoryRef: any;
  modelRef: any;
  model: any;
}
interface PreparedPdbSource {
  label: string;
  dataRef: any;
  trajectoryRef: any;
  modelRef: any;
  structureRef: any;
  componentRef: any;
  representationRef: any;
}
function createStackingConnectorShapeTransformer() {
  return StateTransformer.builderFactory('preparic')({
    name: 'stacking-centroid-connectors',
    display: 'Stacking Centroid Connectors',
    from: SO.Root,
    to: SO.Shape.Provider,
    params: { segments: PD.Value([] as DashSegment[], { isHidden: true }) },
  })({
    apply({ params }) {
      const meshState = MeshBuilder.createState();
      meshState.currentGroup = 0;
      for (const segment of params.segments) {
        addSimpleCylinder(meshState, MolVec3.create(...segment.start), MolVec3.create(...segment.end), {
          radiusTop: 0.045,
          radiusBottom: 0.045,
          radialSegments: 8,
        });
      }
      const mesh = MeshBuilder.getMesh(meshState);
      return new SO.Shape.Provider({
        label: 'Stacking centroid connectors',
        data: mesh,
        params: Mesh.Params,
        getShape: (_ctx, geometry) => Shape.create(
          'Stacking centroid connectors',
          geometry,
          geometry,
          () => Color(0xf4b942),
          () => 1,
          () => 'Centroid-to-centroid stacking interaction',
          undefined,
          1,
        ),
        geometryUtils: Mesh.Utils,
      }, { label: 'Stacking centroid connectors' });
    },
  });
}

let stackingConnectorShapeTransformer: ReturnType<typeof createStackingConnectorShapeTransformer> | undefined;

function getStackingConnectorShapeTransformer() {
  stackingConnectorShapeTransformer ??= createStackingConnectorShapeTransformer();
  return stackingConnectorShapeTransformer;
}

export class MolViewer {
  readonly plugin: PluginContext;
  private modelRef: any = null;
  private model: any = null;
  private structureRef: any = null;
  private loadedCif: string | null = null;
  private loadedId: string | null = null;
  private lastVisibleView: { cif: string; id: string; assemblyId: string } | null = null;
  private activeSource: StructureSource | null = null;
  private stagedSource: StructureSource | null = null;
  private retiredSource: StructureSource | null = null;
  private activePreparedPdbs: PreparedPdbSource[] = [];
  private loadToken = 0;
  private showToken = 0;
  private renderToken = 0;
  private loadQueue: Promise<unknown> = Promise.resolve();
  private showQueue: Promise<unknown> = Promise.resolve();
  private renderQueue: Promise<unknown> = Promise.resolve();
  private clickSubscription: { unsubscribe(): void } | undefined;
  private onPick: (pick: PickInfo) => void = () => undefined;
  private onCifDownloaded: (cif: string) => void = () => undefined;

  private constructor(plugin: PluginContext) {
    this.plugin = plugin;
  }

  static async create(canvas: HTMLCanvasElement, container: HTMLDivElement): Promise<MolViewer> {
    const spec = DefaultPluginSpec();
    const clickBehaviors = new Set([
      PluginBehaviors.Representation.SelectLoci,
      PluginBehaviors.Representation.FocusLoci,
      PluginBehaviors.Camera.FocusLoci,
    ]);
    spec.behaviors = spec.behaviors.filter(behavior => !clickBehaviors.has(behavior.transformer as never));

    const plugin = new PluginContext(spec);
    try {
      await plugin.init();
      const initialized = await plugin.initViewerAsync(canvas, container);
      if (!initialized) throw new Error('WebGL is unavailable in this browser.');
      plugin.canvas3d?.setProps({ renderer: { pickingAlphaThreshold: 0.1 } });

      const viewer = new MolViewer(plugin);
      plugin.managers.interactivity.setProps({ granularity: 'chain' });
      viewer.clickSubscription = plugin.behaviors.interaction.click.subscribe(event => {
        const loci = event.current.loci;

        if (!StructureElement.Loci.is(loci)) return;
        const location = StructureElement.Loci.getFirstLocation(loci);

        if (!location || !Unit.isAtomic(location.unit)) return;
        const labelSeqId = StructureProperties.residue.label_seq_id(location);
        if (!Number.isInteger(labelSeqId) || labelSeqId <= 0) return;
        const pick = {
          labelAsymId: StructureProperties.chain.label_asym_id(location),
          instanceId: location.unit.conformation.operator.instanceId,
          labelSeqId,
        };

        viewer.onPick(pick);
      });
      return viewer;
    } catch (error) {
      plugin.dispose();
      throw error;
    }
  }

  setOnPick(callback: (pick: PickInfo) => void): void {
    this.onPick = callback;
  }

  private downloadedCif: string | null = null;

  setOnCifDownloaded(callback: (cif: string) => void): void {
    this.onCifDownloaded = callback;
    if (this.downloadedCif !== null) {
      try {
        callback(this.downloadedCif);
      } catch {
        // A pending CIF upload must not interrupt rendering.
      }
    }
  }

  setGranularity(granularity: 'chain' | 'residue'): void {
    this.plugin.managers.interactivity.setProps({ granularity });
  }

  loadPdb(id: string, isCurrent: () => boolean = () => true): Promise<AssemblyInfo[] | null> {
    const token = ++this.loadToken;
    const operation = this.loadQueue.then(() => this.loadPdbNow(id, token, isCurrent));
    this.loadQueue = operation.then(() => undefined, () => undefined);
    return operation;
  }

  loadCif(cif: string, id: string, isCurrent: () => boolean = () => true): Promise<AssemblyInfo[] | null> {
    const token = ++this.loadToken;
    this.downloadedCif = null;
    const operation = this.loadQueue.then(() => this.loadCifNow(cif, id, token, undefined, isCurrent));
    this.loadQueue = operation.then(() => undefined, () => undefined);
    return operation;
  }

  clear(): Promise<void> {
    const token = ++this.loadToken;
    this.downloadedCif = null;
    this.showToken += 1;
    this.renderToken += 1;
    const pendingLoad = this.loadQueue;
    const pendingShow = this.showQueue;
    const pendingRender = this.renderQueue;
    const operation = Promise.all([pendingLoad, pendingShow, pendingRender]).then(async () => {
      if (token !== this.loadToken) return;
      await this.plugin.clear();
      this.modelRef = null;
      this.model = null;
      this.structureRef = null;
      this.loadedCif = null;
      this.loadedId = null;
      this.lastVisibleView = null;
      this.activeSource = null;
      this.stagedSource = null;
      this.retiredSource = null;
      this.activePreparedPdbs = [];
    });
    this.loadQueue = operation.then(() => undefined, () => undefined);
    return operation;
  }

  /** Load backend-prepared range PDB text without changing its atom records. */
  showPreparedPdbs(
    ranges: readonly { range: string; pdb: string }[],
    overlays: readonly { label: string; pdb: string }[] = [],
    isCurrent: () => boolean = () => true,
  ): Promise<void> {
    const token = ++this.loadToken;
    this.downloadedCif = null;
    const previousView = this.lastVisibleView;
    const operation = this.loadQueue.then(async () => {
      await this.showQueue;
      await this.renderQueue;
      if (token !== this.loadToken || !isCurrent()) return;
      await this.discardStagedStructure();
      const staged: PreparedPdbSource[] = [];
      let activated = false;
      try {
        for (const item of ranges) {
          if (typeof item.range !== 'string' || !/^range-[1-9]\d*$/.test(item.range)
            || typeof item.pdb !== 'string' || !item.pdb.length) {
            throw new Error('A prepared range PDB is missing its range address or text.');
          }
          const source = await this.stagePreparedPdb(item.pdb, `Prepared ${item.range}`, token, isCurrent);
          if (!source) {
            await this.discardPreparedPdbs(staged);
            return;
          }
          staged.push(source);
        }
        for (const item of overlays) {
          if (typeof item.pdb !== 'string' || !item.pdb.length) continue;
          const source = await this.stagePreparedPdb(item.pdb, item.label, token, isCurrent);
          if (!source) {
            await this.discardPreparedPdbs(staged);
            return;
          }
          staged.push(source);
        }
        if (token !== this.loadToken || !isCurrent()) {
          await this.discardPreparedPdbs(staged);
          return;
        }

        // Keep the current view intact until one guarded state update swaps
        // visibility from the old structures to the fully staged set.
        if (token !== this.loadToken || !isCurrent()) {
          await this.discardPreparedPdbs(staged);
          return;
        }
        const oldRenderRefs = this.renderTaggedRootRefs();
        const oldStructureRef = this.structureRef;
        const oldPreparedPdbs = [...this.activePreparedPdbs];
        const oldActiveSource = this.activeSource;
        const oldRetiredSource = this.retiredSource;
        const state = this.plugin.state.data;
        for (const ref of oldRenderRefs) {
          state.updateCellState(ref, { isHidden: true });
        }
        for (const source of oldPreparedPdbs) {
          if (source.representationRef) {
            state.updateCellState(source.representationRef.ref, { isHidden: true });
          }
        }
        for (const source of staged) {
          state.updateCellState(source.representationRef.ref, { isHidden: false });
        }
        // The readiness check and all visibility updates above run in the
        // same JavaScript turn. No async state commit can apply after a
        // disconnect and hide the last visible view.
        activated = true;
        this.activePreparedPdbs = staged;
        this.modelRef = null;
        this.model = null;
        this.structureRef = null;
        this.loadedCif = null;
        this.loadedId = null;
        this.activeSource = null;
        this.stagedSource = null;
        this.retiredSource = null;
        this.lastVisibleView = previousView;

        // The old refs are hidden now. Cleanup cannot change the visible view
        // if readiness changes while these asynchronous deletions run.
        await this.deleteStateRefs(oldRenderRefs);
        if (oldStructureRef) await this.deleteStateRefs([oldStructureRef]);
        await this.discardPreparedPdbs(oldPreparedPdbs);
        if (oldActiveSource) await this.deleteStructureSource(oldActiveSource);
        if (oldRetiredSource && oldRetiredSource !== oldActiveSource) {
          await this.deleteStructureSource(oldRetiredSource);
        }
        if (staged.length > 0 && token === this.loadToken && isCurrent()) {
          this.plugin.managers.camera.reset();
        }
      } catch (error) {
        if (!activated) await this.discardPreparedPdbs(staged);
        if (token !== this.loadToken || !isCurrent()) return;
        throw error;
      }
    });
    this.loadQueue = operation.then(() => undefined, () => undefined);
    return operation;
  }

  private async stagePreparedPdb(
    pdb: string,
    label: string,
    token: number,
    isCurrent: () => boolean,
  ): Promise<PreparedPdbSource | null> {
    const staged: PreparedPdbSource = {
      label,
      dataRef: null,
      trajectoryRef: null,
      modelRef: null,
      structureRef: null,
      componentRef: null,
      representationRef: null,
    };
    const current = () => token === this.loadToken && isCurrent();
    try {
      staged.dataRef = await this.plugin.builders.data.rawData({ data: pdb, label });
      if (!current()) {
        await this.discardPreparedPdbSource(staged);
        return null;
      }
      staged.trajectoryRef = await this.plugin.builders.structure.parseTrajectory(staged.dataRef as never, 'pdb');
      if (!current()) {
        await this.discardPreparedPdbSource(staged);
        return null;
      }
      staged.modelRef = await this.plugin.builders.structure.createModel(staged.trajectoryRef, { modelIndex: 0 });
      if (!current()) {
        await this.discardPreparedPdbSource(staged);
        return null;
      }
      staged.structureRef = await this.plugin.builders.structure.createStructure(staged.modelRef);
      if (!current()) {
        await this.discardPreparedPdbSource(staged);
        return null;
      }
      const component = await this.plugin.builders.structure.tryCreateComponentStatic(
        staged.structureRef,
        'all',
        { label, tags: [PREPARED_PDB_TAG] },
      );
      if (!component || !current()) {
        await this.discardPreparedPdbSource(staged);
        return null;
      }
      staged.componentRef = component;
      staged.representationRef = await this.plugin.builders.structure.representation.addRepresentation(component, {
        type: 'ball-and-stick',
        color: 'element-symbol',
        typeParams: { visuals: ['intra-bond'] },
      }, { tag: PREPARED_PDB_TAG, initialState: { isHidden: true } });
      if (!current()) {
        await this.discardPreparedPdbSource(staged);
        return null;
      }
      return staged;
    } catch (error) {
      await this.discardPreparedPdbSource(staged);
      if (!current()) return null;
      throw new Error(readableLoadError(error));
    }
  }

  private async discardPreparedPdbSource(source: PreparedPdbSource): Promise<void> {
    await this.deleteStateRefs([
      source.representationRef,
      source.componentRef,
      source.structureRef,
      source.modelRef,
      source.trajectoryRef,
      source.dataRef,
    ]);
  }

  private async discardPreparedPdbs(sources: PreparedPdbSource[]): Promise<void> {
    if (sources.length === 0) return;
    await this.deleteStateRefs(sources.flatMap(source => [
      source.representationRef,
      source.componentRef,
      source.structureRef,
      source.modelRef,
      source.trajectoryRef,
      source.dataRef,
    ]));
  }

  /** Clear stage-4 structures while keeping the selected RCSB assembly restorable. */
  clearPreparedPdbs(): Promise<void> {
    const token = ++this.loadToken;
    const previousView = this.lastVisibleView;
    const pendingLoad = this.loadQueue;
    const pendingShow = this.showQueue;
    const pendingRender = this.renderQueue;
    const operation = Promise.all([pendingLoad, pendingShow, pendingRender]).then(async () => {
      if (token !== this.loadToken) return;
      await this.plugin.clear();
      this.modelRef = null;
      this.model = null;
      this.structureRef = null;
      this.loadedCif = null;
      this.loadedId = null;
      this.activeSource = null;
      this.stagedSource = null;
      this.retiredSource = null;
      this.activePreparedPdbs = [];
      this.lastVisibleView = previousView;
    });
    this.loadQueue = operation.then(() => undefined, () => undefined);
    return operation;
  }

  private async loadPdbNow(id: string, token: number, isCurrent: () => boolean): Promise<AssemblyInfo[] | null> {
    if (token !== this.loadToken || !isCurrent()) return null;
    this.downloadedCif = null;
    await this.showQueue;
    await this.renderQueue;
    if (token !== this.loadToken || !isCurrent()) return null;
    try {
      const data = await this.plugin.builders.data.download({
        url: `https://files.rcsb.org/download/${encodeURIComponent(id.toUpperCase())}.cif`,
        isBinary: false,
        label: id.toUpperCase(),
      });
      if (token !== this.loadToken || !isCurrent()) {
        await this.deleteStateRefs([data]);
        return null;
      }
      if (!data?.obj) throw new Error('The PDB entry could not be downloaded. Check the four-character PDB code.');
      const cif = data.obj.data;
      if (typeof cif !== 'string') throw new Error('The PDB entry did not return a text CIF file.');
      // Keep the last visible structure in place while the network request is
      // pending. Recheck the shared target before replacing it with this CIF.
      if (!cif.trim() || token !== this.loadToken || !isCurrent()) {
        await this.deleteStateRefs([data]);
        return null;
      }
      this.downloadedCif = cif;
      try { this.onCifDownloaded(cif); } catch { /* Saving CIF bytes does not block loading. */ }
      return await this.loadCifNow(cif, id, token, data, isCurrent);
    } catch (error) {
      if (token !== this.loadToken || !isCurrent()) return null;
      throw new Error(readableLoadError(error));
    }
  }

  private async loadCifNow(cif: string, id: string, token: number, downloadedData?: unknown, isCurrent: () => boolean = () => true): Promise<AssemblyInfo[] | null> {
    if (token !== this.loadToken || !isCurrent()) return null;
    if (!cif.trim()) throw new Error('The stored CIF is empty.');
    await this.showQueue;
    await this.renderQueue;
    if (token !== this.loadToken || !isCurrent()) return null;
    await this.discardStagedStructure();
    const staged: StructureSource = {
      cif,
      id: id.toUpperCase(),
      dataRef: null,
      trajectoryRef: null,
      modelRef: null,
      model: null,
    };
    this.stagedSource = staged;
    try {
      const data = downloadedData ?? await this.plugin.builders.data.rawData({ data: cif, label: id.toUpperCase() });
      staged.dataRef = data;
      if (token !== this.loadToken || !isCurrent()) {
        await this.discardStagedStructure(staged);
        return null;
      }
      const trajectory = await this.plugin.builders.structure.parseTrajectory(data as never, 'mmcif');
      staged.trajectoryRef = trajectory;
      if (token !== this.loadToken || !isCurrent()) {
        await this.discardStagedStructure(staged);
        return null;
      }
      const modelRef = await this.plugin.builders.structure.createModel(trajectory, { modelIndex: 0 });
      staged.modelRef = modelRef;
      if (token !== this.loadToken || !isCurrent()) {
        await this.discardStagedStructure(staged);
        return null;
      }
      const model = modelRef.obj?.data;
      if (!model) throw new Error('The shared CIF did not contain a usable model.');
      staged.model = model;
      const assemblies = this.getAssemblies(model);
      this.downloadedCif = cif;
      return assemblies;
    } catch (error) {
      await this.discardStagedStructure(staged);
      if (token !== this.loadToken || !isCurrent()) return null;
      throw new Error(readableLoadError(error));
    }
  }

  activateStagedStructure(): void {
    const staged = this.stagedSource;
    if (!staged?.modelRef || !staged.model) throw new Error('No authorized structure load is ready to display.');
    if (this.retiredSource) void this.deleteStructureSource(this.retiredSource);
    this.retiredSource = this.activeSource;
    this.activeSource = staged;
    this.stagedSource = null;
    this.modelRef = staged.modelRef;
    this.model = staged.model;
    this.loadedCif = staged.cif;
    this.loadedId = staged.id;
  }

  async discardStagedStructure(expected: StructureSource | null = this.stagedSource): Promise<void> {
    if (!expected) return;
    if (this.stagedSource === expected) this.stagedSource = null;
    await this.deleteStructureSource(expected);
  }

  private async deleteStructureSource(source: StructureSource | null): Promise<void> {
    if (!source) return;
    await this.deleteStateRefs([source.modelRef, source.trajectoryRef, source.dataRef]);
  }

  private async deleteStateRefs(refs: unknown[]): Promise<void> {
    const uniqueRefs = [...new Set(refs.filter(ref => ref != null))];
    if (uniqueRefs.length === 0) return;
    try {
      const update = this.plugin.state.data.build();
      for (const ref of uniqueRefs) update.delete(ref as never);
      await update.commit();
    } catch {
      // Cleanup failures do not replace the structure or connectivity error.
    }
  }

  async restoreLastVisibleView(): Promise<void> {
    const snapshot = this.lastVisibleView;
    await this.clearPreparedPdbs();
    if (!snapshot) return;
    const assemblies = await this.loadCif(snapshot.cif, snapshot.id);
    if (!assemblies?.some(assembly => assembly.id === snapshot.assemblyId)) return;
    await this.showStagedAssembly(snapshot.assemblyId, () => true, () => true);
  }

  showAssembly(
    id: string,
    isCurrent: () => boolean = () => true,
    beforeSwap?: () => boolean,
  ): Promise<void> {
    const token = ++this.showToken;
    const operation = this.showQueue.then(() => this.showAssemblyNow(id, token, isCurrent, false, beforeSwap));
    this.showQueue = operation.then(() => undefined, () => undefined);
    return operation;
  }


  showStagedAssembly(
    id: string,
    isCurrent: () => boolean,
    beforeSwap: () => boolean,
  ): Promise<void> {
    const token = ++this.showToken;
    const operation = this.showQueue.then(() => this.showAssemblyNow(id, token, isCurrent, true, beforeSwap));
    this.showQueue = operation.then(() => undefined, () => undefined);
    return operation;
  }

  private async showAssemblyNow(
    id: string,
    token: number,
    isCurrent: () => boolean,
    useStagedSource = false,
    beforeSwap?: () => boolean,
  ): Promise<void> {
    const stagedSource = useStagedSource ? this.stagedSource : null;
    const modelRef = useStagedSource ? stagedSource?.modelRef : this.modelRef;
    const model = useStagedSource ? stagedSource?.model : this.model;
    const current = () => token === this.showToken
      && isCurrent()
      && (!useStagedSource || this.stagedSource === stagedSource);
    if (!current()) return;
    if (!modelRef || !model) throw new Error('Load a PDB entry before choosing an assembly.');
    const validIds = new Set((ModelSymmetry.Provider.get(model)?.assemblies ?? []).map(assembly => assembly.id));
    const isAsymmetricUnit = id === 'model';
    if (!isAsymmetricUnit && !validIds.has(id)) throw new Error(`Assembly ${id} is not available for this entry.`);

    await this.renderQueue;
    if (!current()) return;

    const stagedRefs: unknown[] = [];
    let activated = false;
    try {
      const structure = await this.plugin.builders.structure.createStructure(modelRef, isAsymmetricUnit
        ? { name: 'model', params: {} }
        : { name: 'assembly', params: { id } });
      stagedRefs.push(structure);
      if (!current()) {
        await this.deleteStateRefs(stagedRefs);
        return;
      }

      const component = await this.plugin.builders.structure.tryCreateComponentStatic(
        structure,
        'polymer',
        { label: 'Polymer', tags: [RENDER_TAG] },
      );
      if (component) stagedRefs.push(component);
      if (!current()) {
        await this.deleteStateRefs(stagedRefs);
        return;
      }

      let representationRef: { ref: string } | undefined;
      if (component) {
        representationRef = await this.plugin.builders.structure.representation.addRepresentation(component, {
          type: 'cartoon',
          color: 'chain-id',
        }, { tag: RENDER_TAG, initialState: { isHidden: true } });
        stagedRefs.push(representationRef);
      }
      if (!current()) {
        await this.deleteStateRefs(stagedRefs);
        return;
      }
      if (beforeSwap && !beforeSwap()) {
        await this.deleteStateRefs(stagedRefs);
        return;
      }
      if (useStagedSource) this.activateStagedStructure();

      const oldStructureRef = this.structureRef;
      const retiredSource = this.retiredSource;
      const oldRenderRefs = this.renderTaggedRootRefs().filter(ref => ref !== component?.ref);
      const state = this.plugin.state.data;
      // The target check, share commit, source activation, and visibility
      // change run without an await, so stale work cannot replace the last view.
      for (const ref of oldRenderRefs) state.updateCellState(ref, { isHidden: true });
      if (representationRef) state.updateCellState(representationRef.ref, { isHidden: false });
      this.structureRef = structure;
      activated = true;
      if (this.loadedCif && this.loadedId) {
        this.lastVisibleView = { cif: this.loadedCif, id: this.loadedId, assemblyId: id };
      }
      try { this.plugin.managers.camera.reset(); } catch { /* A camera reset does not undo the committed view. */ }

      // Only hidden refs from the previous assembly are removed after the
      // authorized swap. A stale cleanup cannot blank the new assembly.
      await this.deleteStateRefs(oldRenderRefs);
      if (oldStructureRef) await this.deleteStateRefs([oldStructureRef]);
      if (retiredSource && retiredSource !== this.activeSource) {
        this.retiredSource = null;
        await this.deleteStructureSource(retiredSource);
      }
    } catch (error) {
      if (activated) return;
      await this.deleteStateRefs(stagedRefs);
      if (!current()) return;
      throw new Error(readableLoadError(error));
    }
  }

  listChains(): ChainInfo[] {
    const structure = this.structureRef?.obj?.data;
    if (!structure) return [];

    type RnaAtoms = { o3?: Vec3; phosphorus?: Vec3 };
    const grouped = new Map<string, ChainUnitInput & {
      seqIds: Set<number>;
      rnaBreakAfter: Set<number>;
      residueAtoms: Map<number, RnaAtoms>;
    }>();
    for (const unit of structure.units) {
      if (!Unit.isAtomic(unit) || unit.polymerElements.length === 0) continue;
      const element = unit.elements[0];
      const location = StructureElement.Location.create(structure, unit, element);
      const labelAsymId = StructureProperties.chain.label_asym_id(location);
      const instanceId = unit.conformation.operator.instanceId;
      const key = `${labelAsymId}\u0000${instanceId}`;
      let chain = grouped.get(key);
      if (!chain) {
        chain = {
          labelAsymId,
          instanceId,
          isIdentity: unit.conformation.operator.isIdentity,
          seqIds: new Set<number>(),
          rnaBreakAfter: new Set<number>(),
          residueAtoms: new Map<number, RnaAtoms>(),
        };
        grouped.set(key, chain);
      }
      const hierarchy = unit.model.atomicHierarchy;
      const labelSeqIds = hierarchy.residues.label_seq_id;
      const residueIndices = new Set<ResidueIndex>();
      for (let i = 0; i < unit.polymerElements.length; i += 1) {
        residueIndices.add(unit.residueIndex[unit.polymerElements[i]]);
      }
      for (const residueIndex of residueIndices) {
        const seqId = labelSeqIds.value(residueIndex);
        if (!Number.isInteger(seqId) || seqId <= 0) continue;
        chain.seqIds.add(seqId);
        const atoms = chain.residueAtoms.get(seqId) ?? {};
        const atomPosition = (...names: string[]): Vec3 | undefined => {
          const atomIndex = names.map(name => hierarchy.index.findAtomOnResidue(residueIndex, name))
            .find(index => index >= 0);
          if (atomIndex === undefined) return undefined;
          const point = MolVec3.create(0, 0, 0);
          unit.conformation.position(atomIndex, point);
          return [point[0], point[1], point[2]];
        };
        atoms.o3 ??= atomPosition("O3'", 'O3*');
        atoms.phosphorus ??= atomPosition('P');
        chain.residueAtoms.set(seqId, atoms);
      }
    }
    for (const chain of grouped.values()) {
      const sequence = [...chain.seqIds].sort((a, b) => a - b);
      for (let index = 1; index < sequence.length; index += 1) {
        const before = sequence[index - 1];
        const after = sequence[index];
        if (after !== before + 1) continue;
        const oxygen = chain.residueAtoms.get(before)?.o3;
        const phosphorus = chain.residueAtoms.get(after)?.phosphorus;
        if (!oxygen || !phosphorus) continue;
        const distance = Math.hypot(
          oxygen[0] - phosphorus[0],
          oxygen[1] - phosphorus[1],
          oxygen[2] - phosphorus[2],
        );
        if (!Number.isFinite(distance) || distance > 2.5) chain.rnaBreakAfter.add(before);
      }
    }
    return nameChains([...grouped.values()]);
  }

  render(spec: RenderSpec): Promise<void> {
    const token = ++this.renderToken;
    const operation = this.renderQueue.then(() => this.renderNow(spec, token));
    this.renderQueue = operation.then(() => undefined, () => undefined);
    return operation;
  }

  private async renderNow(spec: RenderSpec, token: number): Promise<void> {
    if (token !== this.renderToken || !this.structureRef) return;
    await this.deleteRenderTaggedNodes();
    if (token !== this.renderToken) return;

    if (!spec.focus) {
      const polymer = await this.plugin.builders.structure.tryCreateComponentStatic(
        this.structureRef,
        'polymer',
        { label: 'Polymer', tags: [RENDER_TAG] },
      );
      if (polymer && token === this.renderToken) {
        await this.plugin.builders.structure.representation.addRepresentation(polymer, {
          type: 'cartoon',
          color: 'chain-id',
        }, { tag: RENDER_TAG });
      }
      if (token === this.renderToken) await this.renderStackingOverlays(spec.stackings ?? [], token);
      return;
    }

    const primaryFocus = spec.range ? rangeExpression(spec.range) : chainExpression(spec.focus);
    const focusExpressions = [primaryFocus, ...(spec.alwaysVisible ?? []).map(rangeExpression)];
    const focus = unionExpression(focusExpressions);
    const focusComponent = await this.plugin.builders.structure.tryCreateComponentFromExpression(
      this.structureRef,
      focus,
      'phase1-focus',
      { label: 'Focused chain or range and pinned ranges', tags: [RENDER_TAG] },
    );
    if (focusComponent && token === this.renderToken) {
      await this.plugin.builders.structure.representation.addRepresentation(focusComponent, {
        type: 'cartoon',
        color: 'chain-id',
      }, { tag: RENDER_TAG });
    }

    const sidechainRanges = spec.sidechainRanges ?? [];
    if (sidechainRanges.length > 0) {
      const sidechainComponent = await this.plugin.builders.structure.tryCreateComponentFromExpression(
        this.structureRef,
        unionExpression(sidechainRanges.map(sidechainRangeExpression)),
        'phase1-sidechain-sticks',
        { label: 'Sidechain sticks for selected ranges', tags: [RENDER_TAG] },
      );
      if (sidechainComponent && token === this.renderToken) {
        await this.plugin.builders.structure.representation.addRepresentation(sidechainComponent, {
          type: 'ball-and-stick',
          color: 'element-symbol',
          typeParams: { visuals: ['intra-bond'] },
        }, { tag: RENDER_TAG });
      }
    }
    if (token !== this.renderToken || spec.others === 'invisible') return;

    const restComponent = await this.plugin.builders.structure.tryCreateComponentFromExpression(
      this.structureRef,
      complementExpression(focus),
      'phase1-rest',
      { label: 'Other chains', tags: [RENDER_TAG] },
    );
    if (!restComponent || token !== this.renderToken) return;
    await this.plugin.builders.structure.representation.addRepresentation(restComponent, {
      type: 'cartoon',
      color: 'chain-id',
      ...(spec.others === 'transparent' ? { typeParams: { alpha: 0.2 } } : {}),
    }, { tag: RENDER_TAG });
  }

  private async renderStackingOverlays(stackings: readonly StackingViewerInput[], token: number): Promise<void> {
    if (stackings.length === 0 || token !== this.renderToken || !this.structureRef) return;

    const participantRanges = uniqueStackingParticipantRanges(stackings);
    if (participantRanges.length > 0) {
      const participants = await this.plugin.builders.structure.tryCreateComponentFromExpression(
        this.structureRef,
        unionExpression(participantRanges.map(rangeExpression)),
        'stage3-stacking-participants',
        { label: 'Selected stacking participants', tags: [RENDER_TAG] },
      );
      if (participants && token === this.renderToken) {
        await this.plugin.builders.structure.representation.addRepresentation(participants, {
          type: 'ball-and-stick',
          color: 'element-symbol',
          typeParams: { visuals: ['intra-bond'] },
        }, { tag: RENDER_TAG });
      }
    }

    if (token !== this.renderToken) return;
    const segments = stackings.flatMap(stacking => stackingDashSegments(
      stacking.rnaCentroid,
      stacking.proteinCentroid,
    ));
    if (segments.length === 0) return;
    await this.plugin.state.data.build()
      .toRoot()
      .apply(getStackingConnectorShapeTransformer(), { segments }, { tags: [RENDER_TAG] })
      .apply(ShapeRepresentation3D, {}, { tags: [RENDER_TAG] })
      .commit();
  }

  private async deleteRenderTaggedNodes(): Promise<void> {
    await this.deleteStateRefs(this.renderTaggedRootRefs());
  }

  private renderTaggedRootRefs(): string[] {
    const tagged = StateSelection.Generators.root.subtree().withTag(RENDER_TAG).select(this.plugin.state.data);
    if (tagged.length === 0) return [];
    const taggedRefs = new Set(tagged.map(cell => cell.transform.ref));
    return tagged
      .filter(cell => !taggedRefs.has(cell.transform.parent))
      .map(cell => cell.transform.ref);
  }

  private getAssemblies(model: any): AssemblyInfo[] {
    const symmetry = ModelSymmetry.Provider.get(model);
    const polymerChains = this.getPolymerChainIds(model);
    const assemblies = symmetry?.assemblies ?? [];
    if (assemblies.length === 0) {
      return [{
        id: 'model',
        label: 'Asymmetric unit',
        participants: [...polymerChains].map(labelAsymId => ({ labelAsymId, count: 1 })),
      }];
    }

    const details = this.getAssemblyDetails(model);
    return assemblies.map(assembly => {
      const counts = new Map<string, number>();
      for (const group of assembly.operatorGroups) {
        const count = group.operators.length;
        for (const asymId of group.asymIds ?? []) {
          if (polymerChains.has(asymId)) counts.set(asymId, (counts.get(asymId) ?? 0) + count);
        }
      }
      return {
        id: assembly.id,
        label: `Assembly ${assembly.id}`,
        details: details.get(assembly.id),
        participants: [...counts].map(([labelAsymId, count]) => ({ labelAsymId, count })),
      };
    });
  }

  private getPolymerChainIds(model: any): Set<string> {
    const ids = new Set<string>();
    const chains = model.atomicHierarchy.chains;
    for (let i = 0; i < chains._rowCount; i += 1) {
      const entityId = chains.label_entity_id.value(i);
      const entityIndex = model.entities.getEntityIndex(entityId);
      const entityType = model.entities.data.type.value(entityIndex);
      if (entityType === 'polymer') ids.add(chains.label_asym_id.value(i));
    }
    return ids;
  }

  private getAssemblyDetails(model: any): Map<string, string> {
    const result = new Map<string, string>();
    if (!MmcifFormat.is(model.sourceData)) return result;
    const table = model.sourceData.data.db.pdbx_struct_assembly;
    for (let i = 0; i < table._rowCount; i += 1) {
      const detail = table.oligomeric_details.value(i)?.trim();
      if (detail) result.set(table.id.value(i), detail);
    }
    return result;
  }

  dispose(): void {
    this.loadToken += 1;
    this.showToken += 1;
    this.renderToken += 1;
    this.clickSubscription?.unsubscribe();
    this.plugin.dispose();
  }
}

function readableLoadError(error: unknown): string {
  const message = error instanceof Error ? error.message : String(error);
  if (/404|not found|no data|empty response/i.test(message)) return 'Entry not found. Check the four-character PDB code.';
  return message || 'The structure could not be loaded. Check your connection and try again.';
}
