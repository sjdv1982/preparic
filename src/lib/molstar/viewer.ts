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
import { chainExpression, complementExpression, rangeExpression, sidechainRangeExpression, unionExpression } from './expressions';
import { stackingDashSegments, type DashSegment } from './stacking-geometry';
import { nameChains, type ChainUnitInput } from '../chains';
import { uniqueStackingParticipantRanges } from '../stackings';
import type { AssemblyInfo, ChainInfo, PickInfo, RenderSpec, StackingViewerInput } from '../types';

const RENDER_TAG = 'preparic-phase1-render';
interface StructureSource {
  cif: string;
  id: string;
  dataRef: any;
  trajectoryRef: any;
  modelRef: any;
  model: any;
}
const StackingConnectorShape = StateTransformer.builderFactory('preparic')({
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
      const cif = data.obj?.data;
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
    await this.discardStagedStructure();
    if (!snapshot) return;
    if (this.retiredSource) {
      const failedSource = this.activeSource;
      if (this.structureRef) {
        await this.deleteStateRefs([this.structureRef]);
        this.structureRef = null;
      }
      this.activeSource = this.retiredSource;
      this.retiredSource = null;
      this.modelRef = this.activeSource.modelRef;
      this.model = this.activeSource.model;
      this.loadedCif = this.activeSource.cif;
      this.loadedId = this.activeSource.id;
      await this.deleteStructureSource(failedSource);
    }
    if (this.structureRef && this.activeSource?.cif === snapshot.cif
      && this.activeSource.id === snapshot.id) return;
    if (!this.activeSource || this.activeSource.cif !== snapshot.cif || this.activeSource.id !== snapshot.id) {
      const assemblies = await this.loadCif(snapshot.cif, snapshot.id);
      if (!assemblies?.some(assembly => assembly.id === snapshot.assemblyId)) return;
      this.activateStagedStructure();
    }
    await this.showAssembly(snapshot.assemblyId);
  }

  showAssembly(id: string): Promise<void> {
    const token = ++this.showToken;
    const operation = this.showQueue.then(() => this.showAssemblyNow(id, token));
    this.showQueue = operation.then(() => undefined, () => undefined);
    return operation;
  }

  private async showAssemblyNow(id: string, token: number): Promise<void> {
    if (token !== this.showToken) return;
    if (!this.modelRef || !this.model) throw new Error('Load a PDB entry before choosing an assembly.');
    const validIds = new Set((ModelSymmetry.Provider.get(this.model)?.assemblies ?? []).map(assembly => assembly.id));
    const isAsymmetricUnit = id === 'model';
    if (!isAsymmetricUnit && !validIds.has(id)) throw new Error(`Assembly ${id} is not available for this entry.`);

    if (this.structureRef) {
      await this.plugin.state.data.build().delete(this.structureRef).commit();
      this.structureRef = null;
    }
    if (token !== this.showToken) return;
    const structure = await this.plugin.builders.structure.createStructure(this.modelRef, isAsymmetricUnit
      ? { name: 'model', params: {} }
      : { name: 'assembly', params: { id } });
    if (token !== this.showToken) {
      await this.plugin.state.data.build().delete(structure).commit();
      return;
    }
    this.structureRef = structure;
    await this.render({ focus: null });
    if (token === this.showToken) {
      this.plugin.managers.camera.reset();
      if (this.loadedCif && this.loadedId) {
        this.lastVisibleView = { cif: this.loadedCif, id: this.loadedId, assemblyId: id };
      }
      const retired = this.retiredSource;
      this.retiredSource = null;
      await this.deleteStructureSource(retired);
    }
  }

  listChains(): ChainInfo[] {
    const structure = this.structureRef?.obj?.data;
    if (!structure) return [];

    const grouped = new Map<string, ChainUnitInput & { seqIds: Set<number> }>();
    for (const unit of structure.units) {
      if (!Unit.isAtomic(unit) || unit.polymerElements.length === 0) continue;
      const element = unit.elements[0];
      const location = StructureElement.Location.create(structure, unit, element);
      const labelAsymId = StructureProperties.chain.label_asym_id(location);
      const instanceId = unit.conformation.operator.instanceId;
      const key = `${labelAsymId}\u0000${instanceId}`;
      let chain = grouped.get(key);
      if (!chain) {
        chain = { labelAsymId, instanceId, isIdentity: unit.conformation.operator.isIdentity, seqIds: new Set<number>() };
        grouped.set(key, chain);
      }
      const labelSeqIds = unit.model.atomicHierarchy.residues.label_seq_id;
      for (let i = 0; i < unit.polymerElements.length; i += 1) {
        const elementIndex = unit.polymerElements[i];
        const residueIndex = unit.residueIndex[elementIndex];
        const seqId = labelSeqIds.value(residueIndex);
        if (Number.isInteger(seqId) && seqId > 0) chain.seqIds.add(seqId);
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
      .apply(StackingConnectorShape, { segments }, { tags: [RENDER_TAG] })
      .apply(ShapeRepresentation3D, {}, { tags: [RENDER_TAG] })
      .commit();
  }

  private async deleteRenderTaggedNodes(): Promise<void> {
    const tagged = StateSelection.Generators.root.subtree().withTag(RENDER_TAG).select(this.plugin.state.data);
    if (tagged.length === 0) return;
    const taggedRefs = new Set(tagged.map(cell => cell.transform.ref));
    const update = this.plugin.state.data.build();
    for (const cell of tagged) {
      if (!taggedRefs.has(cell.transform.parent)) update.delete(cell.transform.ref);
    }
    await update.commit();
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
