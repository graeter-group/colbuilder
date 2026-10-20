# Attachment-preserving ratio replacement

## Configuration

Add these settings to a normal Colbuilder run with the desired geometry and
the original N/C crosslink recipes. The percentage is the fraction REMOVED.

```yaml
replace_bool: true
ratio_replace: 50
ratio_replace_scope: enzymatic
ratio_replace_mode: preserve_attachment
ratio_replace_seed: 14
topology_generator: true
force_field: amber99
```

`ratio_replace_mode` defaults to `random`, preserving the previous seeded
selection algorithm. `preserve_attachment` is an explicit, different sampling
policy. `ratio_replace_seed` is optional; omitting it keeps the existing
time-based seed. The selected seed is always saved in the replacement report.

In constrained mode, `ratio_replace` must be present, including for 0%, and
`manual_replacements` cannot override the selection. Direct PDB replacement
and replacement of an in-memory System use the same selector and contract.
Automatic cleanup of chemically unpaired markers is a separate operation:
an unpaired marker is not a complete crosslink and does not count as attachment.

## Exact guarantee

Let a model be one model/caps ID, not an individual polypeptide chain. For
every model participating in at least one complete crosslink BEFORE ratio
replacement, at least one complete crosslink must remain AFTER replacement.
Initially unlinked models are not made linked by this operation.

Crosslinks are indivisible chemical entities: a PYD comprises LYX + LY2 +
LY3 and two forming bonds, but contributes ONE unit to each incident model's
crosslink degree. A PYD spanning three models protects all three models.

Coverage counts all complete crosslink families, including unselected scopes.
For example, a retained MOLD or glucosepane can keep a model attached when
PYD is removed. This is a model-attachment guarantee, not equivalence of
trivalent and divalent chemistry. To restrict which families can provide
attachment would require another explicitly defined policy.

The per-scope removal targets retain the existing rule:

`round(number_of_eligible_entities_in_scope * ratio_replace / 100)`.

Rounding therefore still applies when a percentage is not an integer count;
Python's nearest-integer/ties-to-even behavior is unchanged. N and C quotas
are separate when the configured terminal recipes identify `enzymatic_n`
and `enzymatic_c`. Unclassified enzymatic entities remain one `enzymatic`
scope; the selector cannot infer missing terminal metadata.

All scopes are solved jointly with binary variables (one per removable
entity), using `scipy.optimize.milp`/HiGHS. Quotas are equalities. For each
initially linked model, removals cannot exceed its initial degree minus one.
A seeded random objective chooses between feasible solutions; this is NOT
uniform sampling of all feasible solutions. No global Python RNG state is
modified. Reproduction assumes the same input and solver/software versions.

Infeasible quotas raise an error before ratio mutations are sent to Chimera.
The percentage is never silently lowered. The error also reports the maximum
total number removable under the attachment constraints without scope quotas;
that diagnostic upper limit does not assert feasibility of a particular N/C
split or percentage after rounding.

## PYD50 example

For 20 PYD arranged on ten disjoint model pairs, each with one N and one C
site, remove at most one PYD per pair. At 50%, five N removals and five C
removals must affect complementary sets of model pairs. Ten complete PYD
remain and all twenty initially linked models retain one crosslink.

The previous unconstrained PYD50 run removed both sites on 14--70, 23--58
and 56--63. Correct marker counts and intact backbone bonds did not prevent
these new losses of attachment. This mode specifically prevents that outcome.

## Validation and reports

1. The selector records the complete removed and retained identities, seed,
   target quotas and model degrees before/after.
2. After Chimera, the usual complete-marker and backbone checks still run.
   The remaining entity identities and model coverage must additionally match
   the constrained plan. The contract survives network refresh/regrouping.
3. Before AMBER exports its final TOP/GRO, an independent pass rereads every
   final ITP. Original model/residue/chain identities are mapped through the
   registered source-residue ordering, not guessed from repeated residue IDs.
4. All inter-model bonds must match the planned forming bonds exactly, once
   each. A PYD counts as retained only if both forming bonds exist. Missing,
   duplicated, unexpected or split crosslinks, missing protected models and
   mapping inconsistencies stop topology generation.

The existing validations of angles, proper/improper dihedrals, 1-4 pairs,
marker templates, backbone continuity and GRO atom ordering remain active.

Outputs in both the working data and final topology directory:

- `replacement_report.json`: seed, policy, requested/realized integer quotas,
  complete removed/retained identities, initial/final entity degree per model.
- `crosslink_network.json`: final chemical graph plus the replacement report.
- `replacement_attachment_validation.json`: confirmation derived from actual
  final ITP bonds, protected IDs and their observed degrees. Produced only
  after the constrained AMBER validation passes.

A topology-only run needs the original replacement report as well as its
matching geometry. A bare post-replacement PDB cannot reveal which models
were linked before removal; requesting a guarantee without that provenance
fails rather than silently reconstructing a weaker condition. The AMBER
builder accepts the in-memory contract or `replacement_report.json` in its
working directory, and the public pipeline copies it from the chosen caps
directory when constrained mode is requested.

The final ITP validator currently supports AMBER99 only. A constrained run
requesting Martini topology is rejected explicitly rather than exporting an
unchecked topology. Geometry-only replacement can still produce the PDB and
its contract for subsequent supported topology generation.

## What this does not guarantee

- Attachment of every individual chain in a triple helix.
- Retention of both terminal sites on every model pair.
- Preservation of every original edge, connected component or load path in
  a general graph: deleting a bridge can split a component without producing
  a zero-degree model. The requested policy prevents zero degree only.
- Absence of sliding, unfolding, group detachment or instability under load.

Consequently this is not a certification of mechanical stability, and the
constrained ensemble is scientifically distinct from unconstrained random
crosslink depletion. Document the policy when comparing simulations.

## Tests

`tests/test_replacement_policy.py` covers 100 seeded PYD50 selections, exact
quotas, infeasibility, preservation of legacy random draws, RNG isolation,
unselected AGE attachment, small general graphs checked by exhaustive search,
configuration errors, contract persistence and deliberately corrupted ITPs.

The native regression script exercises real Chimera replacement, backbone
coordinate preservation, AMBER/pdb2gmx, final ITP validation and `grompp`:

```bash
python tests/run_ratio_topology_regression.py \
  --case /absolute/path/to/saved_pyd_case \
  --output /absolute/path/to/new_validation_directory \
  --mode preserve_attachment --seed 14
```

It copies the saved pre-replacement caps; it does not overwrite the source
case. No minimization, equilibration or pulling is performed by this test.

Local validation on 2026-09-07: 415 tests passed, including cleanup before
constrained selection. Both the saved-caps regression and a complete CLI
geometry/replacement/topology run passed for Rattus PYD50 with seed 14:
5 N + 5 C PYD, 45 ITPs, all 20 initially linked models still attached,
and successful `grompp` without `-maxwarn`. The full CLI run also preserved
139421 N/CA/C/O coordinates across the 30 marker-to-LYS mutations.
