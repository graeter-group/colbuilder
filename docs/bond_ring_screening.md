# Bond/ring screening and local repair

## Why this exists

The saved MOLD+glucosepane template contained SER518.B CB-OG passing through
the five-membered AGS523.A ring before mixing. Correct covalent terms and
finite forces did not identify this bad initial geometry. Repeating the
template produced ten copies of the defect. A subsequent screen of *all*
heavy bonds also found ten AGS494.B/HYP265.C crossings in that old fibril.

The new code does not add exclusions, stiffen bonds, change atom types, or
rewrite the force field. It distinguishes conformer construction from MD.

## Geometric definition

`core/utils/ring_geometry.py` builds an ordered minimum cycle basis from
each residue's RTP heavy-atom bond graph, using NetworkX. AGS has the expected
five- and seven-membered fused rings; standard cyclic residues and the other
crosslink markers are screened as well. Global cycles spanning multiple
residues/crosslinks are not interpreted as chemical ring surfaces.

For each ring, nearby heavy covalent bond segments are screened against a
centroid triangle fan. A second test checks the polygon projected onto the
best-fit plane. Bonds sharing a ring atom are excluded because their contact
with the surface is intrinsic to the ring, not penetration by another bond.

A confirmed penetration must satisfy both tests, have one unique surface
intersection, and be separated from the ring perimeter and segment endpoints
by more than 0.02 A. Endpoints must be on opposite sides of the plane by
more than that tolerance. The report includes out-of-plane deformation,
endpoint signed distances, intersection position and involved atoms.

Grazing, coplanar and inconsistent fan/plane contacts are classified as
ambiguous. Missing ring atoms and unsupported residue templates are reported,
not silently certified as clear. This is a nonperiodic geometric screen, not
a proof of topological linking or of dynamical stability. Strongly puckered
rings have no unique spanning surface; ambiguous cases require inspection.
PBC handling remains a separate simulation-setup concern.

## Sequence optimization

`core/sequence/steric_optimization.py` uses the same ring screen on both full
translated helices. Its spatial indices cache fixed bonds/rings; each trial
only reevaluates surfaces affected by the moved atoms.

- Confirmed penetrations add a large objective penalty and disqualify a
  candidate even when forming bonds and nonbonded distances would pass.
- Ambiguous contacts add a softer penalty and remain explicit warnings.
- The existing differential-evolution search changes only permitted acyclic
  sidechain torsions. No bond within an AGS ring is rotated or stretched.
- When expanding the flexible neighborhood, residues on penetrating bonds
  are considered before the worst atom-pair contact.
- Replay and the final check after all optimization jobs use the same rule.
- `optimize_crosslinks.py` also writes `iteration_N.rings.json`, next to the
  optimization intermediate PDBs, for a full post-optimization two-copy scan.

An exhausted search for a valid current-site conformer follows the existing
optimization-error/retry mechanism. Ordinary compressed-contact warnings
and their previous thresholds have not been made stricter. Unrelated,
unchanged defects in the full post-optimization scan are reported explicitly.

## Mixing

`core/geometry/ring_validation.py` reads all selected caps as one coordinate
system. It builds native residue bonds, peptide bonds respecting the existing
chain/TER boundaries, and the explicit crosslink bonds from the chemical
network when available. Therefore a ring and the crossing bond may belong
to different models or future ITP groups.

In `shared_sites`, screening takes place after grafting, before publishing
the staged caps. Local repair is attempted only for rotatable standard
sidechains on confirmed crossing bonds. It never moves a crosslink marker,
any backbone atom, a proline/hydroxyproline ring, or an entire model.

Single-torsion residues such as SER use a deterministic angular scan; more
complex sidechains use a seeded differential-evolution search. Candidates
are scored against the entire fixed surrounding structure, with exclusions
derived from bonds rather than distances. A move must reduce confirmed
crossings, introduce no new confirmed crossing, and not worsen an already
catastrophic nearest heavy-atom contact. Ring ambiguity is also penalized.
The complete structure is checked again after PDB coordinate rounding.

This intentionally limited repair cannot solve every defect. For example,
penetration by a fixed HYP ring is not repaired by distorting that ring.
Such cases require regenerating the AGE conformer or a separately designed
relaxation procedure, and remain unresolved in the report.

Caps publication remains transactional in `shared_sites`: if strict validation
fails, no composed caps or model-state changes are published. The diagnostic
report is retained. Backbone records, TER boundaries and atom serials remain
unchanged by local repair. Marker coordinates and forming bonds remain fixed.

The legacy `whole_models` strategy performs a final screen but does not apply
native-sidechain repair, preserving whole-model coordinate selection semantics.

## Configuration and reports

```yaml
mix_strategy: shared_sites
mix_ring_repair: true
ring_penetration_policy: warn
```

- `mix_ring_repair`: defaults to true; attempts the restricted local repair
  after shared-site composition. False performs detection only.
- `ring_penetration_policy`: defaults to `warn`, preserving the user-selected
  warning-only export policy. `error` rejects unresolved confirmed crossings
  after mixing. It does not turn ordinary distance warnings into errors and
  is not a general certification of MD-ready geometry.

`ring_geometry_report.json` is written in the mixing working directory and
included under `ring_geometry` in `crosslink_mix_report.json`. It records
before/after crossings, atom identities, attempted successful repairs,
sidechain displacement and incomplete screening information.

`core/topology/coordinate_validation.py` repeats the heavy-bond/ring screen
on the actual GRO and ITP bonds after topology generation. This catches
changes due to export/rounding and includes bonds across residue boundaries.
It writes `<gro_basename>.rings.json` and warns without changing coordinates
or the existing topology-stage warning-only policy.

## Tests and operational workflow

`tests/test_ring_geometry.py` covers the real AGS523/SER518 fixture,
cross-model contacts, graph-derived fused rings, confirmed/absent/grazing
contacts, moving-ring and moving-bond paths, strict no-write failure,
nonfinite coordinates, and preservation of marker/backbone coordinates and
native bond lengths during repair. Existing mixing and sequence tests remain
part of the full regression suite.

`tests/run_ring_caps_regression.py` copies a saved caps directory into a new
output directory before running the repair and invariance checks. With
`--policy error`, unresolved non-repairable defects are expected to fail;
`--policy warn` retains their diagnosis and allows inspecting partial repair.

For a new MOLD+GCP system, start from MOLD without AGS/LGX, rerun the AGE
addition YAMLs in sequence, and then rerun mixing. Do not feed the previously
pathological AGE residues back as an allegedly repaired template. Validate
the final GRO/ITPs, prepare the intended simulation box, and independently
check EM/equilibration. Passing this screen is not a substitute for those steps.
