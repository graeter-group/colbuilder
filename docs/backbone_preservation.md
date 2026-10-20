# Backbone preservation during crosslink replacement

For the step-by-step implementation walkthrough in Italian, see
[backbone_preservation_it.md](backbone_preservation_it.md).

The subsequent fixes for ratio selection and final crosslink grouping are
documented separately in [crosslink_network.md](crosslink_network.md).

## Failure addressed

Chimera can interpret an elongated peptide C-N contact as a chain break and write
a new `TER` record when a crosslink residue is converted to LYS or ARG. For
example, replacing LYX103 with LYS in the saved rat PYD 50% replacement case
introduced five breaks before GLY104, despite unchanged C/N coordinates.

With its normal chain separation rules, `pdb2gmx` then creates CLYS/NGLY termini
instead of a peptide bond. This changes atoms, charges, bonded terms and load
transmission. `grompp` can still compile that chemically unintended system.

## Behaviour and configuration

Backbone preservation is automatic. No new YAML key is required. Existing
`ratio_replace`, `ratio_replace_scope`, `replace_file`, `manual_replacements`,
and `contact_distance` values retain their meanings.

The protection applies to system-based replacement, direct-PDB replacement,
manual/auto-fix replacement and the mixer's Chimera mutation step. It protects
both LYS and ARG replacements, irrespective of the selected crosslink scope.
It does not change the ratio selection algorithm or fix other mixing issues.

The original PDB is the authority for polymer identity:

- Original residue order, chain IDs, insertion codes, TER boundaries and caps
  are retained. A genuine TER within the same chain ID is not removed.
- All non-target residues are copied from the original, not from Chimera.
- Target backbone atom records and coordinates are copied from the original,
  with only the residue name changed. Missing backbone atoms in Chimera output
  are therefore restored, including the previously observed missing oxygen.
- Only the requested new sidechain atoms come from Chimera. The expected LYS
  or ARG heavy-atom inventory and finite coordinates are checked.
- Atom and TER serials are renumbered. Untouched CONECT edges and explicit
  backbone edges are remapped; CONECT edges involving replaced sidechains are
  dropped. The force-field residue templates define the new sidechain bonds.
- Alternative conformers, ambiguous/repeated residue identities, multiple
  models in a mutation input, missing targets, conflicting mutations, missing
  residues, reordered identities, incomplete sidechains or an incomplete input
  target backbone cause an explicit failure. These inputs are not guessed at.
- Backbone displacement above 0.01 A in a mutated residue's Chimera output is
  rejected to avoid grafting a sidechain from a changed coordinate frame.
  Smaller formatting differences do not propagate: original coordinates are
  always retained.

No distance-based rule is used to decide whether to restore a peptide link.
An already strained peptide bond remains covalently connected if the original
PDB places its residues in the same polymer segment. This change does not move
atoms closer together, modify force constants or replace minimization.

## Transaction and topology validation

`core/geometry/backbone.py::preserve_backbone` snapshots the complete affected
PDB files before invoking Chimera. It validates every output before publishing
the reconstructed files, and restores every affected input if the mutation,
validation or writing fails. Files are replaced atomically, one file at a time;
this is not a multi-file transaction against power loss or an OS process kill.

After AMBER ITP generation, `validate_backbone_topology` compares the ordered
PDB residues and their polymer segments to the ITP. It requires the expected
N-CA, CA-C, C-O and inter-residue C-N bonds. Legitimate terminal oxygen renaming
is allowed only at input segment ends. Peptide bonds across real input
boundaries are rejected. A failure is raised as `TopologyGenerationError`
(`TOP_ERR_005`), so the group cannot silently be skipped and exported as an
apparently successful, smaller system.

The PDB preservation also protects inputs to Martini generation, but the new
bond-by-bond ITP validator is specific to all-atom AMBER, not coarse-grained
topology validation.

The legacy `snapshot_residue_atoms` / `repair_missing_backbone_atoms` helpers
remain available for compatibility. They are not sufficient to preserve
connectivity and are no longer used by the replacement pipeline.

## Existing broken outputs

Restart replacement from the **pre-replacement** geometry. If an input already
contains an artificial TER from an earlier run, the guard cannot distinguish
it from an intentional break and must preserve it.

Regenerate the PDB, GRO, ITP/TOP and position-restraint files together. Do not
add a C-N bond to an existing CLYS/NGLY topology: that would leave incorrect
terminal atoms, charges and interaction terms in place. Do not globally strip
TER records or impose `-chainsep id`, because genuine breaks may share a chain ID.

This fix leaves the initial geometry unchanged. Long crosslink bonds, distorted
HYP carbonyls, steric clashes, minimization, equilibration and pulling stability
are separate checks.

## Tests

Run the unit/regression suite with the project's development dependencies:

```bash
python -m pytest -q tests
```

`tests/test_backbone_preservation.py` tests artificial versus genuine TERs,
chain boundaries, unchanged coordinates, missing backbone oxygen, insertion
codes, serial/CONECT handling, input ambiguity, transactional rollback and
rejection of missing or spurious peptide bonds in AMBER output.

An opt-in end-to-end reproducer replays the saved replacement instructions of
the rat PYD case using real Chimera and GROMACS, preserves the saved group
partition, and rebuilds/validates every group with the normal chain separation:

```bash
conda run --no-capture-output -n colbuilder python \
  tests/run_backbone_replacement_regression.py \
  --case /absolute/path/to/pyd_50pp \
  --output /absolute/path/to/a/new/regression_directory
```

The case needs `.tmp/geometry_gen/T`, the saved replacement instructions, and
the previous topology output. The output directory must not exist. Original
campaign data are never overwritten. `regression_result.json` and GROMACS logs
record the results; the zero-step TPR is diagnostic, not a production MD input.
