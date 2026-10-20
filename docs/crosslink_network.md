# Replacement and final crosslink topology: shared chemical identities

## Scope

This change fixes duplicate ratio sampling and inconsistent final molecular
grouping. It retains the separate backbone preservation guard documented in
[backbone_preservation.md](backbone_preservation.md). No new YAML key is needed.

The prior backbone change did not modify the selector or lattice cutoff, but
its original end-to-end reproducer replayed a saved mutation list and group
partition. That test missed defects in selecting and grouping a fresh run.
The new regressions exercise those stages too.

## Defects reproduced

The connect file represented the same pair as both `11 18` and `18 11`.
The selector sorted each row but did not deduplicate chemical entities before
sampling. Twenty PYD crosslinks became forty candidates. In the observed run,
sampling half of that duplicated list removed sixteen distinct PYD, generating
48 mutations instead of 30. Deduplicating instructions after sampling cannot
correct this bias.

Separately, `model.connect` used a 3 A lattice-growth contact test, while AMBER
searched for forming atoms within 5 A only after grouping. Consequently,
partners separated by 3.25 A could be exported in separate molecule types and
never reach bond generation. Some named PYD bonds were longer than 5 A as well.
Changing the growth cutoff globally would also alter the generated fibril.

## One chemical inventory

`core/geometry/crosslink_network.py` defines:

- `Marker`: model id, residue id/insertion code, chain, residue name and atoms.
- `Entity`: one complete crosslink, its scope, constituent residues and exact
  forming-atom bonds. PYD contains three residues and two inter-marker bonds.
- `Network`: the complete inventory used for selection, grouping and AMBER.

The supported families have atom-specific rules. In particular:

- PYD: LYX:C12--LY2:CB and LYX:C13--LY3:CG;
- glucosepane: LGX:CE--AGS:NZ;
- MOLD: LZD:CE--LZS:NZ1.

Different divalent families cannot be matched to each other merely because
their atoms are nearby. Marker reuse across entities is rejected.

### Assignment and long bonds

Generic pairing requires unique compatible forming atoms within 5 A, on
different models. The same bound determines these final topology candidates;
the lattice-growth constant remains 3 A and is not reused for final grouping.

For a trivalent terminal crosslink, a configured N/C positional recipe can
identify a complete trio even with a longer forming bond, but only when all
three prescribed residues are unambiguous in the central marker's local caps
neighborhood. The neighborhood uses the existing 5 A marker-atom contact
criterion. This is a chemical-identity rule, not permission to join arbitrary
nearest residues or distant models. Multiple assignments or unpaired markers
cause an explicit error rather than an arbitrary choice.

Long recipe-assigned bonds retain their normal force-field parameters and are
reported with their distances. Coordinates are not moved. These warnings mean
that the initial geometry needs relaxation; they are not a claim of MD safety.
Without recipe provenance, a bond beyond the normal cutoff is not inferred.

## Replacement

`replace_in_system` now resolves entities from only the active caps. Stale
model files are not counted. `select_replacements` samples unique complete
entities, separating N, C and non-enzymatic buckets according to the selected
scope. The existing nearest-integer rule is retained for non-integral targets.
The reported random seed makes the selection reproducible.

For 10 N-terminal and 10 C-terminal PYD, a 50% ratio means exactly 5 removals
per terminal: 10 entities and 30 residue mutations. AGE selection uses complete
pairs; LGX becomes LYS and AGS becomes ARG. An AGE-only ratio leaves PYD intact.
Within a mixed AGE bucket the ratio is over the selected scope, not a promise
to remove the same percentage of every chemical subtype independently.

Manual instructions remain authoritative, but removal of only part of a known
complete entity is rejected before Chimera. Instructions targeting already
unpaired markers can still remove those markers. Chain case is resolved
consistently with the backbone guard. After mutation, the code checks the exact
expected residual marker inventory and rejects unexpected changes or orphans.

`preserve_backbone` still surrounds Chimera. The new chemical checks do not
strip TER records globally, alter backbone coordinates or change native bond
parameters. Direct-PDB ratio replacement also uses complete chemical entities.

The legacy connect-based helper remains for external callers; it now
deduplicates normalized groups and entity identities and rejects overlapping
assignments rather than counting the same marker twice.

## Final molecular groups

For an optional model-attachment constraint on ratio replacement, see
[Attachment-preserving replacement](attachment_preserving_replacement.md).
The default `random` selector is unchanged; `preserve_attachment` solves
scope quotas jointly and checks the final AMBER ITP bonds against the plan.

`install_network` derives connected model groups from the surviving entities,
updates `System.model.connect`, reloads marker objects from the final caps
without applying lattice transforms a second time, and writes the same graph
to disk. It does not overwrite the pre-replacement connect file.

`build_amber99` refreshes this validated inventory against the actual caps.
Without a prior inventory, it resolves one from the current caps and config.
The grouped models and the list of bonds therefore come from the same graph.
Mixed A/B caps are copied and merged according to each model's ownership, not
the first model's directory. Caps outside the active System are ignored.

The resulting molecule types are graph components of models; they need not
be a single covalent component at atom level. A collagen triple helix still
has three polypeptide chains unless covalent links connect them.

## Exact mapping and fail-closed export

During merging AMBER records the original model/residue identity and its
position in the merged residue order. The ITP atom index is resolved by this
identity plus the exact forming-atom name, rather than selecting the closest
C12/C13 by coordinates. This avoids ambiguity and GRO serial wrapping in large
merged groups. The coordinate fallback is retained only for legacy callers
that did not use the registered merge path.

Bond completion errors are propagated as `TopologyGenerationError` instead of
being logged and ignored. `core/topology/crosslink_validation.py` independently
checks the final ITP for:

- exactly the expected inter-marker bonds, with no duplicates or wrong roles;
- all forming atoms and exactly one correct partner for each forming site;
- every angle and proper dihedral spanning a crosslink bond;
- all required 1-4 pairs at shortest graph distance three, following the
  existing AMBER hydrogen-pair convention;
- the PYD LY2:CB--LYX:C12--LYX:C13--LY3:CG path and its endpoint 1-4 pair;
- marker RTP atom types, charges, internal bonds and internal impropers;
- absence of additional explicit exclusions in these crosslinked AMBER ITPs.

The backbone validator remains a separate check. Missing group files, failed
merges and missing GRO atoms cannot silently reduce the exported system.
GRO atom counts and atom-name order are checked against their ITPs on assembly.
GROMACS still resolves bonded parameters; the regression scripts additionally
run `grompp` without `-maxwarn` to check compilation with the actual force field.

## Reports

- `.tmp/replace_crosslinks/replacement_report.json`: seed, ratio, scope and
  available/removed/remaining entity counts in each bucket.
- `crosslink_network.json`: exact surviving marker identities, bonds and their
  initial distances. Written in the replacement/topology work directories and
  copied to the final topology directory.

These files document the graph used for the run; an unrelated stale connect
file cannot change the graph after validation.

## Regressions

Run the automated suite with `PYTHONPATH=src python -m pytest -q tests` from
the repository root, using an interpreter with pytest and the project
dependencies installed. On the validation machine, the test interpreter was
`/opt/anaconda3/bin/python`; real external-tool runs used the `colbuilder`
Conda environment and its editable installation of this repository.

- `test_crosslink_network.py`: 100 seeds for PYD50, percentage boundaries,
  complete-entity removal, reversed/duplicate connect rows, invalid ratios,
  long recipe-defined PYD, ambiguous/missing partners, AGE chemistry and scope,
  unchanged growth cutoff and consistent in-memory/on-disk groups.
- `test_crosslink_validation.py`: missing/duplicate/wrong bonds, angles,
  proper dihedrals and 1-4 pairs, exact mapping, failed completion and export.
- Existing backbone and Martini tests remain in the suite.

Real Chimera/GROMACS reproducers, to run in the colbuilder environment with
new output directories:

```bash
python tests/run_ratio_topology_regression.py --case /path/to/pyd_50pp --output /new/output
python tests/run_saved_system_regression.py --case /path/to/10_GCP --ratio 30 --output /new/output
python tests/run_saved_system_regression.py --case /path/to/pyd_50pp_MLD_50pp --output /new/output
python tests/run_full_pipeline_regression.py --config /path/to/pyd_50pp/config.yaml --output /new/output
```

The first script runs actual selection, Chimera mutation, final grouping and
the public topology pipeline, rather than replaying old instructions/groups.
The last script also reruns sequence and geometry generation from the YAML.
These tests do not establish equilibrated mechanics, force-field scientific
accuracy or universal support for every input. The strict final ITP checks are
AMBER-specific; full Martini dynamics and the entire mixing workflow are not
certified by these all-atom regressions.
