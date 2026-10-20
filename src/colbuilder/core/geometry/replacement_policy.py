"""Joint ratio selection with an optional minimum model-attachment constraint."""

from collections import Counter

import numpy as np

from .crosslink_network import CrosslinkIntegrityError


def model_incidence(entities):
    """Count complete chemical entities, not their individual forming bonds."""
    return Counter(mid for entity in entities
                   for mid in {port[0] for bond in entity.bonds for port in bond})


def entity_record(entity):
    return {"kind": entity.kind, "scope": entity.scope, "markers": entity.markers}


def record_identity(record):
    return record["kind"], tuple(sorted(tuple(key) for key in record["markers"]))


def constrained_selection(network, buckets, targets, rng):
    """Solve binary removal quotas and model coverage jointly, without rounding down."""
    from scipy.optimize import Bounds, LinearConstraint, milp
    from scipy.sparse import csc_matrix

    candidates = [entity for key in sorted(buckets)
                  for entity in sorted(buckets[key], key=lambda e: e.identity)]
    if not candidates:
        return []
    scopes = sorted(buckets)
    degree = model_incidence(network.entities)
    models = sorted(degree)
    # x_i = 1 removes an entire entity. Unselected scopes remain in the degree
    # budget, so an existing AGE can also maintain a model's attachment.
    rows, cols = [], []
    scope_row = {key: i for i, key in enumerate(scopes)}
    model_row = {mid: len(scopes) + i for i, mid in enumerate(models)}
    for j, entity in enumerate(candidates):
        rows.append(scope_row[entity.scope])
        cols.append(j)
        for mid in model_incidence([entity]):
            rows.append(model_row[mid])
            cols.append(j)
    matrix = csc_matrix((np.ones(len(rows)), (rows, cols)),
                        shape=(len(scopes) + len(models), len(candidates)))
    quotas = [targets[key]["removed"] for key in scopes]
    lower = np.array(quotas + [0] * len(models), dtype=float)
    upper = np.array(quotas + [degree[mid] - 1 for mid in models], dtype=float)
    bounds = Bounds(np.zeros(len(candidates)), np.ones(len(candidates)))
    integral = np.ones(len(candidates))
    result = milp(np.array([rng.random() for _ in candidates]), integrality=integral,
                  bounds=bounds, constraints=LinearConstraint(matrix, lower, upper),
                  options={"mip_rel_gap": 0.0})
    if result.status == 2:
        maximum = milp(-np.ones(len(candidates)), integrality=integral, bounds=bounds,
                       constraints=LinearConstraint(matrix[len(scopes):],
                                                    lower[len(scopes):], upper[len(scopes):]),
                       options={"mip_rel_gap": 0.0})
        limit = (f" At most {int(round(-maximum.fun))}/{len(candidates)} eligible "
                 "crosslinks can be removed even without the N/C quota constraints."
                 if maximum.status == 0 else "")
        raise CrosslinkIntegrityError(
            "Infeasible preserve_attachment replacement: requested removals "
            f"{dict(zip(scopes, quotas))} would leave initially linked models "
            "without crosslinks. Quotas have not been reduced; no mutations were selected."
            + limit
        )
    if result.status != 0 or result.x is None:
        raise CrosslinkIntegrityError(f"Replacement constraint solver failed: {result.message}")
    chosen = np.rint(result.x).astype(int)
    totals = matrix @ chosen
    if (np.any(np.abs(result.x - chosen) > 1e-6)
            or np.any(totals < lower) or np.any(totals > upper)):
        raise CrosslinkIntegrityError("Replacement solver returned an invalid discrete selection")
    return [entity for entity, remove in zip(candidates, chosen) if remove]


def attachment_report(before, retained):
    initial, final = model_incidence(before), model_incidence(retained)
    return {
        "protected_models": sorted(initial),
        "before_degree": {str(mid): initial[mid] for mid in sorted(initial)},
        "after_degree": {str(mid): final[mid] for mid in sorted(initial)},
        "newly_unlinked_models": sorted(initial.keys() - final.keys()),
    }


def validate_attachment_plan(report, entities):
    """Keep the replacement contract through geometry and topology generation."""
    if not report or report.get("mode") != "preserve_attachment":
        return
    degree = model_incidence(entities)
    missing = sorted(mid for mid in report["attachment"]["protected_models"] if not degree[mid])
    if missing:
        raise CrosslinkIntegrityError(
            f"Replacement attachment contract violated: models {missing} lost all crosslinks"
        )
    planned = Counter(record_identity(record) for record in report["retained_crosslinks"])
    actual = Counter(entity.identity for entity in entities)
    if actual != planned:
        raise CrosslinkIntegrityError("Final crosslink entities differ from the replacement plan")


def validate_replacement_mode(config):
    """Do not silently bypass the constrained mode through manual precedence."""
    if getattr(config, "ratio_replace_mode", "random") != "preserve_attachment":
        return
    if getattr(config, "ratio_replace", None) is None:
        raise CrosslinkIntegrityError("preserve_attachment requires ratio_replace, including 0")
    if getattr(config, "manual_replacements", None):
        raise CrosslinkIntegrityError(
            "preserve_attachment is a ratio-selection mode and cannot be combined with manual_replacements"
        )
