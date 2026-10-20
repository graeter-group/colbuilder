# Shared-site crosslink mixing

## Two intentionally separate strategies

- `whole_models` (default): the existing model-component assignment. A selected
  model comes entirely from one template; additional sites need not survive
  selecting a different template. Existing ratio assignment is unchanged.
- `shared_sites` (opt-in): mix alternative chemistries at coincident physical
  residue loci and preserve the union of complete additional crosslinks.

`ratio_mix: "A:50 B:50"` associates A/B with the first/second `files_mix` entry.
In shared-sites mode it means 50% of each *alternative crosslink family*, not
50% of all models, all crosslinks, marker residues, or atoms. A zero-weight
variant can still donate additional sites: their retention is independent of
the alternative-site ratio.

## Matching and counting

1. Generate, cap, and clean boundary orphans for each homogeneous variant.
2. Resolve complete chemical entities using the existing crosslink network
   resolver. Incomplete or ambiguous entities are errors, not a mixing choice.
3. Match by model ID, chain and residue number (including insertion code).
   Identical entities are common and retained once. Different chemistries
   sharing at least two residues, with one complete footprint contained in the
   other, are alternatives. Thus PYD's three-marker footprint can match MOLD's
   two-marker subset. The released third PYD residue comes from the native
   residue already present in the MOLD template.
4. Disjoint, additional entities are retained once from their donor variant.
   One-residue conflicts, partial overlap or multiple competing matches fail
   clearly; they are not silently interpreted as compatible extra crosslinks.
5. Group alternatives by their complete positional recipes, excluding lattice
   model IDs. N and C recipes therefore receive separate quotas. This does
   not depend on enzymatic/non-enzymatic classification: MOLD is still an AGE.
6. Convert percentages into integer per-family quotas by largest remainder;
   ties follow `ratio_mix` label order. For 10 N and 10 C sites, 50:50 means
   5 A + 5 B at each terminal. With three sites, 50:50 necessarily becomes
   2 A + 1 B, and both target and achieved counts appear in the report.
7. Select whole alternative-site backbone components to meet every quota
   exactly. Additional GCP connections do not constrain this selection.
   If the quotas cannot be met without splitting a backbone component, fail
   explicitly instead of silently changing the ratio. The search is
   deterministic by model ID and does not consume global random state;
   it is not a random ensemble of spatial arrangements.

## Geometry and backbone integrity

Each model uses the complete backbone of its selected variant. Models not
in alternative loci use the first template. Extra marker sidechains are
transferred from the donor template with **one proper rigid transform per
complete crosslink**, fitted to N/CA/C across all its marker residues.
Separate marker transforms would deform the crosslink, so are not used.

`mix_alignment_tolerance` (default 0.25 A, allowed >0 to <=1 A) is the maximum
N/CA/C fit residual. This is a compatibility check, not a bond-length target
or a steric-clash cutoff. An incompatible donor fails rather than pulling
backbones together. It should not be increased merely to force an output.

All backbone coordinates, residue order and TER/MODEL boundaries remain
exactly those of the chosen base model. Non-crosslink sidechains are retained.
Source marker sidechains are copied as a rigid entity. Input cap files are
read-only; composed caps are staged and validated before publication under
`.tmp/mixing_crosslinks/_mixed_sites/`. Existing composed caps are not
overwritten. Use a fresh working directory for a new independent run.

Validation rejects changed native sequence, different clipped model sets or
backbone residue boundaries, missing/extra markers, changed partner identity,
or forming-bond length changes above PDB rounding tolerance (0.004 A).
Final connectivity is rebuilt from the complete composed entity ledger,
including additional GCP edges between differently assigned models. The
existing topology builder consumes this ledger and completes bonded terms
and 1-4 pairs; shared-site mixing does not introduce force-field parameters
or distance-based nonbonded exclusions.

This does not certify MD stability. Sidechain transfer can create new contacts
with surrounding residues, and source geometry may already contain clashes.
Inspect the existing coordinate warnings, then minimize and equilibrate.

## Reports and related bug fixes

`crosslink_mix_report.json` records requested ratios, positional recipes,
integer quotas, achieved counts, selected donor for every entity, model
backbone provenance and rigid-fit residuals. It is written in the mixing
directory and copied to the exported topology folder.
`crosslink_network.json` records the final entities, bonds and distances.

Two fixes apply to both mixing workflows:

- Orphan instructions already contain native targets. `ARG` must remain ARG,
  not pass through marker-to-native conversion a second time and become LYS.
  The existing backbone-preserving Chimera wrapper remains in use.
- Final whole-model orphan detection reads exactly the selected caps, rather
  than deduplicating same-named files across A/B directories. Any mutations
  are sent to the selected file's actual directory. Shared-site composition
  instead validates complete retained entities and never silently deletes
  additional sites after composition.

## Tests

`tests/test_shared_site_mixing.py` covers matching, per-family ratios and
rounding, common/additional retention, ratios 0/100, component coupling and
infeasible quotas, rigid geometry and handedness, exact backbone preservation,
failure rollback, configuration, active-file ownership and native ARG swaps.
Existing backbone, replacement and topology tests must also remain passing.

For a completed CLI run with `debug: true`, the integration check is:

```sh
python tests/run_shared_site_mixing_regression.py \
  --case /path/to/completed_shared_sites_run \
  --output /path/to/new_validation_directory \
  --evaluate-forces
```

This compares every composed backbone atom against its selected source,
validates each ITP against its merged PDB, runs `gmx grompp` without `maxwarn`,
and optionally evaluates the initial forces in a diagnostic box padded by 2 nm
(the exported GRO bounding box is not a production solvent box). A zero-step evaluation is not
minimization or an MD stability test. The `--independent-audit` option accepts
an external audit module exposing `audit_case(label, case_path)` for an
additional graph-completeness audit.
