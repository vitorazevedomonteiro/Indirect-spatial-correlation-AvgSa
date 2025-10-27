"""
utils package
-------------
This package contains modular utility functions used for:
1. GMM prediction generation
2. Standard deviation and correlation calculations
3. Final correlation synthesis

Modules:
    - gmm_calculator.py
    - stdev_and_corr.py
    - stdev_and_corr_total.py
    - correlation_final.py
    - correlation_final_total.py
"""

from .gmm_calculator import compute_gmm_predictions
from .stdev_and_corr import (
    process_stdev_combinations,
    compute_correlations,
    compute_numerators,
    compute_denominators,
)
from .stdev_and_corr_total import (
    process_stdev_combinations_total,
    compute_numerators_total,
    compute_denominators_total,
)
from .corr_final import compute_final_corr
from .corr_final_total import compute_final_corr_total

__all__ = [
    "compute_gmm_predictions",
    "process_stdev_combinations",
    "compute_correlations",
    "compute_numerators",
    "compute_denominators",
    "compute_final_corr",
    "process_stdev_combinations_total",
    "compute_correlations_total",
    "compute_numerators_total",
    "compute_denominators_total",
    "compute_final_corr_total"
]