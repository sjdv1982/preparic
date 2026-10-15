import { StructureElement, StructureProperties, Unit } from 'molstar/lib/mol-model/structure.js';
import { MmcifFormat } from 'molstar/lib/mol-model-formats/structure/mmcif.js';
import { ModelSymmetry } from 'molstar/lib/mol-model-formats/structure/property/symmetry.js';
import { PluginBehaviors } from 'molstar/lib/mol-plugin/behavior.js';
import { PluginContext } from 'molstar/lib/mol-plugin/context.js';
import { DefaultPluginSpec } from 'molstar/lib/mol-plugin/spec.js';
import { StateSelection } from 'molstar/lib/mol-state/index.js';
import { chainExpression, complementExpression, rangeExpression, unionExpression } from './expressions';
import { nameChains, type ChainUnitInput } from '../chains';
import type { AssemblyInfo, ChainInfo, PickInfo, RenderSpec } from '../types';

const RENDER_TAG = 'preparic-phase1-render';

export class MolViewer {
  readonly plugin: PluginContext;
  private modelRef: any = null;
  private model: any = null;
  private structureRef: any = null;
  private loadToken = 0;
  private showToken = 0;
  private renderToken = 0;
  private loadQueue: Promise<unknown> = Promise.resolve();
  private showQueue: Promise<unknown> = Promise.resolve();
  private renderQueue: Promise<unknown> = Promise.resolve();
  private clickSubscription: { unsubscribe(): void } | undefined;
  private onPick: (pick: PickInfo) => void = () => undefined;

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

  setGranularity(granularity: 'chain' | 'residue'): void {
    this.plugin.managers.interactivity.setProps({ granularity });
  }

  loadPdb(id: string): Promise<AssemblyInfo[] | null> {
    const token = ++this.loadToken;
    const operation = this.loadQueue.then(() => this.loadPdbNow(id, token));
    this.loadQueue = operation.then(() => undefined, () => undefined);
    return operation;
  }

  clear(): Promise<void> {
    const token = ++this.loadToken;
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
    });
    this.loadQueue = operation.then(() => undefined, () => undefined);
    return operation;
  }

  private async loadPdbNow(id: string, token: number): Promise<AssemblyInfo[] | null> {
    if (token !== this.loadToken) return null;
    await this.showQueue;
    await this.renderQueue;
    if (token !== this.loadToken) return null;
    await this.plugin.clear();
    this.modelRef = null;
    this.model = null;
    this.structureRef = null;
    try {
      const data = await this.plugin.builders.data.download({
        url: `https://models.rcsb.org/${encodeURIComponent(id.toLowerCase())}.bcif`,
        isBinary: true,
        label: id.toUpperCase(),
      });
      if (token !== this.loadToken) return null;
      const trajectory = await this.plugin.builders.structure.parseTrajectory(data, 'mmcif');
      if (token !== this.loadToken) return null;
      const modelRef = await this.plugin.builders.structure.createModel(trajectory, { modelIndex: 0 });
      if (token !== this.loadToken) return null;
      const model = modelRef.obj?.data;
      if (!model) throw new Error('The PDB entry did not contain a usable model.');
      this.modelRef = modelRef;
      this.model = model;
      return this.getAssemblies(this.model);
    } catch (error) {
      if (token !== this.loadToken) return null;
      this.modelRef = null;
      this.model = null;
      throw new Error(readableLoadError(error));
    }
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
    if (token === this.showToken) this.plugin.managers.camera.reset();
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
      return;
    }

    const focusExpressions = [chainExpression(spec.focus), ...(spec.alwaysVisible ?? []).map(rangeExpression)];
    const focus = unionExpression(focusExpressions);
    const focusComponent = await this.plugin.builders.structure.tryCreateComponentFromExpression(
      this.structureRef,
      focus,
      'phase1-focus',
      { label: 'Focused chains and pinned ranges', tags: [RENDER_TAG] },
    );
    if (focusComponent && token === this.renderToken) {
      await this.plugin.builders.structure.representation.addRepresentation(focusComponent, {
        type: 'cartoon',
        color: 'chain-id',
      }, { tag: RENDER_TAG });
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
