# Export gens fixtures

`1b7f.snapshot.json` is a normalized, pure planner input derived from §10 of
`plans/export-gens-plan.md`: clip 2/1; the four selected anchors at nucleotides
4, 5, 6 and 11; and the range-2/range-3 matrix over fragments 3–10. Its PDB
byte members are synthetic placeholders because the gens planner does not read
or render them. `$bytes_base64` tags are decoded by the fixture loader before
the snapshot is passed to the planner.

`1b7f.catalogue.golden.json` is the expected §5.2 catalogue projection derived
from the written R0–R12 recipes. It records every physical gen path, mode,
range and direction, action-directory inventory, final pool, and external
pool contract with its possible catalogue producers. It is independent of the
planner implementation and does not read the user's current mounts.
