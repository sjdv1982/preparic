# Web app first phase outline

The overall function of the web app is to define selections on a molecular structure. The current outline is for the first phase: there will be follow-up phase where there is interaction with a server. The app has stages: the follow-up phases will add additional stages.

Throughout the stages, there will be a Mol* viewer with all buttons disabled: only camera and selection is working.

The first stage is to select a PDB code and a biological unit. As soon as a valid PDB code is entered in a text field, Mol* will load it and a drop-down menu will appear with the biological units (index + participating chains) . By default, the first biounit will be selected. The selected biounit is visualized in the viewer. When the user presses OK, the next stage happens. In general, the user can always go back a stage but this erases the data of the current stage.

The second stage is to select ranges within the chains of the biounit. A range consists of (chain, start, end). An Excel-style table is shown for the ranges: the user can add and remove rows, move them up or down in the table, and click on a cell to have it selectable in Mol*: the chain cell makes click-in-viewer select a chain, the start and end cell make it select a residue, but only residues within chain in the chain cell. Note that in the current row, start can be larger than end: in that case, as soon as the user selects something else than a cell in the current row, start and end are swapped.

Note that "chain" is duplicated if the chain occurs multiple times in the biounit: chain A becomes A (for the identity matrix), A-1, A-2, ... . Note that Mol* has some troubles with this, see ./molstar-bio-unit.md.

When a chain cell is being selected, the molecule is rendered normally (cartoon for all chains). If a start/end cell is being selected, only the chain in question is being rendered. A drop-down button selects what to do with the other chains: visible (i.e. normal), transparent, or invisible. In addition, every range row has a checkbox next to it, if ticked that makes the range row always visible.  
