# AWSEM gamma tables

`AWSEM_2015.json` and `AWSEM_membrane_2015.json` are the published AWSEM-MD
(2015) contact and burial parameters, taken from the `frustratometer` Python
package (Jaafari et al.), `frustratometer/data/`.

Format, which `glueforge/potential.py` parses:

- `gamma_array` — 1260 floats, laid out in the order of `segment_definition`
- `segment_definition` — `Burial (3, 20)`, `Direct (1, 20, 20)`,
  `Protein (1, 20, 20)`, `Water (1, 20, 20)`; 60 + 400 + 400 + 400 = 1260
- `alphabet` — `A R N D C Q E G H I L K M F P S T W Y V`, the AWSEM/LAMMPS
  ordering, **not** alphabetical. The loader reorders onto this package's
  alphabetical convention.

Not implemented here: the electrostatic term, which AWSEM applies separately
with its own coupling constant and a 10 A screening length. The contact and
burial terms are the ones the frustration index is built on.
