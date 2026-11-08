# Preparic stages

This describes what the six stages do, what each stage's Continue requires,
and which files hold the project. The requirements and their checks are in
Alaric's `plans/final/oracle-contract.md`. The stage model is in
`plans/final/state-transitions.md`, and the complete transformer and
shared-cell inventory is in `plans/final/contract-boundaries.md` (§2 and §8).
Identifiers such as T5 or S11 below refer to that inventory.

## How the stage is decided

All tabs share one project and one stage. The stage is derived only from the
stored decisions, never from whether a computation has finished:

| Stage | Reached when this decision is active and matches its upstream decisions |
|---|---|
| 1 | always |
| 2 | ranges (`ranges.yaml`) |
| 3 | stacking selection (`stacking-selection.json`) |
| 4 | threshold table (`best-fit-threshold-edits.json`) |
| 5 | interaction matrix (`interaction-matrix.json`) |
| 6 | Gaussian classes (`gaussian-classes.json`) |

Each of these decisions carries a basis (the PDB code and the decisions it
was made for) and an `active` flag. **Continue** checks the stage's gate
again immediately before it writes, writes any prerequisite decision first,
and then activates the next stage's decision. **Back** deactivates the
current stage's decision and keeps it. A later Continue restores it when its
basis still applies; otherwise a fresh default replaces it. Computation never
moves the stage, and Continue never waits for a computation to finish before
it writes.

The clip (`best-fit-clip.json`) and the energy-threshold table
(`energy-threshold-edits.json`) are stage parameters, not stage decisions:
they have a basis but no active flag.

A result is fresh when the inputs it echoes equal the current decisions,
including edits not yet acknowledged by the backend. Stale results are shown
greyed or as pending. A gate that needs a result stays closed until it is
fresh. Errors are shown in the stage they belong to.

## Stage 1: structure

Enter a PDB code and choose a biological assembly. The CIF is downloaded first.
Only after a successful download does the tab store the code and the CIF
bytes. The backend reads the entry id back from the stored CIF (T26), and
Continue 1→2 requires it to equal the stored code. A failed download changes
neither. In stages 1–3 the viewer shows the RCSB structure.

## Stage 2: ranges

Define residue ranges. Range 1 is the RNA; the others are protein ranges.
Ranges may be picked on any copy of the assembly: the backend applies the
assembly's operators, so the range is extracted, prepared and shown where it
was picked. An RNA range that contains an unobserved residue is refused, in
the browser and in the backend, with the missing residues named.

The backend extracts the ranges, prepares them with aareduce (RNA:
`--rna --heavy --nalib --no-renumbering`; protein: heavy-atom PDB2PQR with the
original numbering) and detects stackings.

**Continue 2→3** needs confirmed valid rows and a fresh, nonempty stacking
table.

## Stage 3: stacking selection

Choose the stacking interactions that will anchor the RNA. The default
selects every stacking with distance **≤ 5 Å** and slip angle φ **< 60°**,
with no tilt condition. Choices kept from an earlier visit are merged in by
stacking key.

On entering stage 3 the backend starts the best fit of every dinucleotide to
the fragment library (A1), the reference fragments (A2), the pair and base
tables (A3) and their rounded thresholds (A4). It also lists the fragments
(T24) and copies the all-atom range PDBs for the later viewer (T25).

**Continue 3→4** needs a selection keyed exactly by the fresh stacking table,
with at least one stacking selected, a fresh fragment list and fresh all-atom
range PDBs. The fit need not have finished.

## Stage 4: best fits, pair thresholds and the clip

From stage 4 on, the viewer shows only the backend's all-atom range PDBs,
byte for byte as aareduce wrote them, with the reference and best-fit
fragments drawn on top. It never shows the RCSB structure or an old PDB; it
shows a pending message until fresh PDBs arrive. The nucleotide shown in the
stage-4 viewer and sidebar is shared between tabs (S26), but it is not
stored and resets on every transition into or out of stage 4.

**Pair thresholds.** The editable table holds one row per adjacent fragment
pair (`down`, `up`, trinucleotide), with an ovRMSD and a cRMSD threshold. These
become the `grow` thresholds of the export. The table opens empty. The
computed values (A4: the raw value rounded to the smallest multiple of 0.01
strictly above it) are shown beside it. **Copy** fills the table from them,
only while they are fresh, and asks before it overwrites typed values. Values
are finite, positive, with at most two decimals; extra decimals are rejected,
not rounded. Values above 10 Å get a warning. Each cell can be reverted, and a
count shows how many differ from the computed values.

**Fragment information.** For every fragment the stage shows its residues,
sequence, best-fit conformer and fit RMSD, read only.

**Base thresholds are not editable.** The `anchor-refe` thresholds are the
computed per-nucleotide base ovRMSD values of A4, used as they are. A stored
threshold file that still contains `bases` is malformed and loads as
undecided.

**The clip.** `clip1` removes fragments from the 5′ end and `clip2` from the 3′
end. Clipped fragments, and every pair that touches one, are greyed. Their
thresholds stay editable but are not required. Clipping changes only what is
selected, displayed and handed on: the fit, the tables and the gridless
energies are always computed on the whole RNA, and a clip edit never reruns
them. Gaussian jobs, scores and classes may recompute for the kept fragments.

A selected stacking on nucleotide N has up to two **candidate fragments**:
fragment N−1 (where N is the second nucleotide) and fragment N (where N is the
first). **The clip must keep at least one candidate of every selected
stacking.** For example, for a stacking on nucleotide 2, fragments 1 and 2 are
candidates: a 5′ clip of 1 is valid, and the stacking continues through
fragment 2. A 5′ clip of 2 removes both candidates and is invalid until that
stacking is deselected in stage 3. The 3′ end is symmetric. A stacking on a
terminal nucleotide has only one candidate, which must stay. With several
selected stackings, each one limits the clip.

The clip is edited only in stage 4, with debounced writes, and is shown read
only in stages 5 and 6. Continue 3→4 resets a clip whose PDB/RNA-range basis
is stale. It keeps a saved clip that leaves one candidate per selected
stacking. Only when a newly selected stacking has no kept candidate does it
reduce each clip end to its limit. It waits for that write before it opens the
table. This automatic, possibly silent clamp is the only automatic correction
of a decision.

**Continue 4→5** needs every kept pair's thresholds filled and valid
(uncovering a clipped blank pair closes the gate again), a valid clip, the
threshold validator (T10) not failed, fresh computed thresholds (S9) and fresh
gridless energetics.

## Stage 5: energetics and the interaction matrix

The backend scores every reference and best-fit fragment against each protein
range separately, gridless (A5), with one scorer per range. Changing one
protein range rescores only that range. Energies are shown for every fragment
and range, coloured up to about ±20, including clipped rows. A scoring error
on a fragment is shown on its row.

**The interaction matrix** has one row per kept fragment and one column per
protein range. True cells mark the fragments that interact with that range.
Each selected stacking gives an anchor column on its range, and its **locks**
are its candidate fragments that are kept by the clip. Every lock must be
true. An anchor survives the clipping of one of its candidates; its other
candidate remains its lock.

**Recompute** starts from the locks and grows each run while the reference or
best-fit energy has an absolute value of at least 3. It stops before a
low-energy, error or clipped row, and before another same-range anchor's
lock. Runs on one range that touch or overlap are joined; runs on different
ranges may overlap. "Until next anchor" is offered only for the next
downstream anchor on the same range, and is forced on when their locks touch
or overlap. The matrix is read back and validated on every load: shape,
extent, locks, an anchor in every true run, and no true cell on a row with a
scoring error. A failed read disables the editor, names the failed check and
offers Recompute. Nothing is repaired automatically.

**A scoring error on a lock** refuses every matrix for that anchor. The user
either deselects the stacking in stage 3, or clips the failing candidate in
stage 4 while another candidate stays. A lock is never dropped automatically.

**Continue 4→5 reuses a saved matrix** when the selected stackings are
unchanged, adapting it once to the current clip. Clipped rows and their
locks are dropped, and newly uncovered locks are restored as true cells, even
at low energy. Other newly uncovered rows are grown only at a run end that
touched the old clip boundary, with Recompute's rule. Interior edits and ends
the user had stopped short are kept. Otherwise Continue starts from
Recompute. A clip edit that races Continue makes the result invalid and the
project stays in stage 4.

**Energy thresholds** form a second editable table, one value per fragment and
protein range. It opens empty, with the computed values shown beside it: the
smallest multiple of 0.01 strictly above the best-fit energy plus 0.05 (A6
rounding plus 0.05), so each lies more than 0.05 and at most 0.06 above that
energy. The margin absorbs the effect of 0.001 Å coordinate precision, which
can move a best fit on a steep repulsive contact by several thousandths.
Copy fills it while those are fresh, asking before it overwrites.
Values have either sign and at most two decimals. The table is editable in
stage 5 only; stage 6 may show its saved and computed values read-only.
Changing one protein range resets only that range's values; a change of PDB
code or RNA range empties the table. Editing it reruns only its validator
(T27); editing the matrix reruns only the interaction hand-off (T15). Stage 5
does not show the S11 pair-threshold table.

**Continue 5→6** needs the matrix to pass read validation against fresh
energetics, a fresh fit and the current clip; fresh best-fit stacking
geometry; fresh Gaussian scores whose golden checks passed; and a valid
energy threshold for every true matrix cell, with T27 not failed.

## Stage 6: Gaussian stacking filter

The backend scores every kept best-fit candidate and every reference stacking
with the Gaussian stacking model (T28), from Alaric's
`gaussian_stacking_filter.json`. It checks the model's golden vectors on every
run. The browser does no Gaussian scoring of its own. Acceptance uses the
JSON's mask (distance < 5 Å and φ < 60°, strict) and the five-decimal cutoffs
with `score >= cutoff`.

Each selected stacking gets one computed class from its **worst kept
candidate**: the smallest of 50, 60, 70, 80, 90 and 95 whose cutoff it reaches.
A stacking below RRM-95 is shown red and keeps class 95. For classes 90 and 95,
the tilt30 and disp45 prefilters are proposed when every kept candidate has
θ < 30° or φ < 45°, respectively. The reference score is shown for
information only. A clipped candidate contributes nothing.

The class table holds one editable row per selected stacking, with Copy,
per-row revert and difference marks, and writes debounced by 500 ms. It is the
only editable decision in stage 6. S11 pair thresholds may appear as a
read-only comparison, and S24 energy thresholds may appear read-only; neither
has inputs or Copy/Revert controls there. A model or fit change updates the
computed values and the marks, never the user's classes. Stage 4 remains the
only stage where S11 can be edited.

**Validation.** The backend runs Alaric's all-filter check (A9, T21) once for
every kept candidate against all twelve filters. Class edits do not rerun
it; they only re-merge its results (T22). Each row then shows green, amber,
red or blocked from its worst kept candidate. A result implied by a broader
filter is shown as implied, never as 0 %. A failing clipped candidate no
longer affects its stacking; uncovering it restores its contribution.

## Export

Once stage 6 and its upstream results are current, Preparic offers **Download
export ZIP**. The backend validates one snapshot of the decisions and computed
results, then builds the DATA files, manifest and generator catalogue. A
pending or failed producer, an invalid decision, or a hard admission block
keeps the export unavailable. The state graph explains pending and failed
producers; `export_admission` reports addressed decision blocks and warnings.

The ZIP contains one `alaric_project/` folder. Its `DATA/` has the exact RNA
reference and prepared range files, schema-1 constraints and energy
thresholds, kept reference/best-fit fragments, and a manifest with the
snapshot digest and generator inventory. The other folders catalogue
alternative anchor, grow, merge and common-conformer gens for fixed and
unfixed runs. They are not all one runnable project: select compatible gens,
then move the chosen actions from their nested `gen-*` folders into the
project root alongside `DATA/`. The catalogue itself is ignored by Alaric's
project loader until those action directories are moved. Choose only one
conflicting alternative for a pool name; the exporter does not choose merge
locations or remove run-final dead ends. There is no bundled README.

Some valid snapshots carry warnings when the nearest-to-native structure does
not pass a filter. In that case, Preparic shows the reasons and the exact
**“The nearest-to-native structure doesn't pass the filters. Proceed anyway”**
checkbox. Acknowledgment applies only to that snapshot digest. The download
request carries no consent flag, and the browser releases the ZIP only when
the manifest digest still matches the click-time snapshot. Later edits mark a
previous download stale.

Gaussian `anchor` actions retain each selected stacking's current filter.
`anchor-refe` actions use A4's computed base threshold for each unclipped
candidate, deduplicated within its generator; an over-limit rounded value is
written as `ovrmsd: -1` with a comment, as a deliberately disabled positive
control. Kept edited pair thresholds populate `constraints.json`; the edited
energy thresholds of true kept matrix cells populate
`energy-thresholds.json`. Fixed-mode scoring applies the required range
filters in a tandem AND chain; each unfixed range-rooted chain scores only
against its own range.

## Files in `deployment-mounts/`

| File | Content | Kind |
|---|---|---|
| `pdbcode.txt` | PDB code (S1) | decision |
| `biounit.json` | biological assembly (S2) | decision |
| `cif.cif` | downloaded structure (S3) | decision |
| `ranges.yaml` | ranges (S4) | stage-2 decision |
| `stacking-selection.json` | selected stackings (S8) | stage-3 decision |
| `best-fit-threshold-edits.json` | edited pair thresholds (S11) | stage-4 decision |
| `best-fit-clip.json` | `clip1`, `clip2` (S12) | stage-4 parameter |
| `interaction-matrix.json` | interaction matrix (S14) | stage-5 decision |
| `energy-threshold-edits.json` | edited energy thresholds (S24) | stage-5 parameter |
| `gaussian-classes.json` | edited Gaussian classes (S18) | stage-6 decision |
| `stackings.json` | detected stackings (S7) | computed output |
| `best-fit-thresholds.json` | computed RMSD thresholds with raw values (S9) | computed output |
| `computed-energy-thresholds.json` | computed energy thresholds (S15) | computed output |
| `interactions.json` | interaction hand-off (T15) | computed output |
| `gaussian-validation.json` | Gaussian validation (S19) | computed output |
| `pdbs/`, `pdbs-aa/` | extracted and aareduced range PDBs | computed output |
| `clipped-fragments/` | kept reference/best-fit PDBs, kept pair thresholds, `clip.json` (T11) | computed output |

Decisions are read on startup and written by the browser. Computed outputs
are written by the backend and never read back.
