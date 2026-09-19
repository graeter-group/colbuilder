"""
Constants for the ColBuilder Pipeline

This module defines constants used by crosslink optimization in the ColBuilder pipeline,
providing a centralized and consistent reference to ensure uniform behavior across the system.

Key Features:
--------------
1. **Optimization Constants**:
   - `MAX_OPTIMIZATION_ATTEMPTS`: Maximum number of optimization attempts for crosslinks.
   - `MAX_TRIVALENT_DISTANCE`: Maximum allowable distance for trivalent crosslinks (in Ångstroms).
   - `MAX_DIVALENT_DISTANCE`: Maximum allowable distance for divalent crosslinks (in Ångstroms).
   - `CRITICAL_DISTANCE_THRESHOLD`: Distance threshold beyond which optimization is considered critical.

Usage:
------
This module is intended to be imported wherever these constants are required to ensure consistency
and avoid hardcoding values.

Example:
--------
```python
from colbuilder.core.utils.constants import MAX_OPTIMIZATION_ATTEMPTS

if attempts > MAX_OPTIMIZATION_ATTEMPTS:
    raise ValueError("Exceeded maximum optimization attempts.")
```
"""

from typing import Final

# Optimization constants
MAX_OPTIMIZATION_ATTEMPTS: Final = 3
MAX_TRIVALENT_DISTANCE: Final = 7.0
MAX_DIVALENT_DISTANCE: Final = 5.0
CRITICAL_DISTANCE_THRESHOLD: Final = 9.0
