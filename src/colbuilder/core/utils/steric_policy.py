"""Geometric overlap thresholds, not force-field contact distances (values in A).

The sequence optimizer rejects almost-coincident heavy atoms. The final GRO
screen reports finite overlaps as warnings only; export does not certify MD stability.
"""

CATASTROPHIC_HEAVY_A = 0.25
CATASTROPHIC_HYDROGEN_A = 0.15
WARN_HEAVY_A = 2.0
WARN_14_A = 1.5
WARN_HYDROGEN_A = 0.6
