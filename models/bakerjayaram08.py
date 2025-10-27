"""
Created by Vitor Monteiro, 23/09/2025

Compute correlation coefficients for Sa(T), for
NGA-W1 database for periods between [0.01-10]s.

For more details please see:
Baker, J. W., & Jayaram, N. (2008). Correlation of
spectral acceleration values from NGA ground motion models.
Earthquake Spectra, 24(1), 299–317. https://doi.org/10.1193/1.2857544
"""

import numpy as np


def corrBJ08(T1, T2):

    T_min, T_max = min(T1, T2), max(T1, T2)

    C1 = 1 - np.cos((np.pi / 2) - 0.366 *
                    np.log(T_max / np.maximum(T_min, 0.109)))

    if T_max < 0.2:
        C2 = 1 - 0.105 * (1 - (1 / (1 + np.exp(100 * T_max - 5)))
                          ) * ((T_max - T_min) / T_max - 0.0099)
    else:
        C2 = 0

    C3 = C2 if T_max < 0.109 else C1

    C4 = C1 + 0.5 * (np.sqrt(C3) - C3) * (1 + np.cos((np.pi * T_min) / 0.109))

    # Determine the rho based on the conditions
    if T_max < 0.109:
        rho = C2
    elif T_min > 0.109:
        rho = C1
    elif T_max < 0.2:
        rho = np.min([C2, C4])
    else:
        rho = C4

    return rho
