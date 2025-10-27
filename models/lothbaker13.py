"""
Created by Vitor Monteiro, 23/09/2025

Compute cross-spatial correlation coefficients between Sa(T) 
for NGA-W1 database database for periods between [0.01-10.0]s.

For more details please see: 
Loth, C., & Baker, J. W. (2013). 
A spatial cross‐correlation model of spectral accelerations at multiple periods. 
Earthquake Engineering & Structural Dynamics, 42(3), 397–417. https://doi.org/10.1002/eqe.2212
"""


import numpy as np

def CrossSpatialCorrLB13(T1, T2, h):
    # Verify validity of input arguments
    if min(T1, T2) < 0.01:
        raise ValueError('The periods must be greater or equal to 0.01s')
    if max(T1, T2) > 10:
        raise ValueError('The periods must be less or equal to 10s')
    if h < 0:
        raise ValueError('The separation distance must be positive')

    # Tlist
    Tlist = [0.01, 0.1, 0.2, 0.5, 1, 2, 5, 7.5, 10.0]

    # Special case: perfect correlation only if periods are equal and distance is zero
    if T1 == T2 and h == 0:
        return 1.0

    # Model coefficients
    B1 = np.array([
        [0.29, 0.25, 0.23, 0.23, 0.18, 0.1, 0.06, 0.06, 0.06],
        [0.25, 0.30, 0.2, 0.16, 0.1, 0.04, 0.03, 0.04, 0.05],
        [0.23, 0.20, 0.27, 0.18, 0.1, 0.03, 0.0, 0.01, 0.02],
        [0.23, 0.16, 0.18, 0.31, 0.22, 0.14, 0.08, 0.07, 0.07],
        [0.18, 0.10, 0.1, 0.22, 0.33, 0.24, 0.16, 0.13, 0.12],
        [0.10, 0.04, 0.03, 0.14, 0.24, 0.33, 0.26, 0.21, 0.19],
        [0.06, 0.03, 0.0, 0.08, 0.16, 0.26, 0.37, 0.3, 0.26],
        [0.06, 0.04, 0.01, 0.07, 0.13, 0.21, 0.3, 0.28, 0.24],
        [0.06, 0.05, 0.02, 0.07, 0.12, 0.19, 0.26, 0.24, 0.23]
    ])

    B2 = np.array([
        [0.47, 0.4, 0.43, 0.35, 0.27, 0.15, 0.13, 0.09, 0.12],
        [0.4, 0.42, 0.37, 0.25, 0.15, 0.03, 0.04, 0.0, 0.03],
        [0.43, 0.37, 0.45, 0.36, 0.26, 0.15, 0.09, 0.05, 0.08],
        [0.35, 0.25, 0.36, 0.42, 0.37, 0.29, 0.2, 0.16, 0.16],
        [0.27, 0.15, 0.26, 0.37, 0.48, 0.41, 0.26, 0.21, 0.21],
        [0.15, 0.03, 0.15, 0.29, 0.41, 0.55, 0.37, 0.33, 0.32],
        [0.13, 0.04, 0.09, 0.2, 0.26, 0.37, 0.51, 0.49, 0.49],
        [0.09, 0.0, 0.05, 0.16, 0.21, 0.33, 0.49, 0.62, 0.6],
        [0.12, 0.03, 0.08, 0.16, 0.21, 0.32, 0.49, 0.6, 0.68]
    ])

    B3 = np.array([
        [0.24, 0.22, 0.21, 0.09, -0.02, 0.01, 0.03, 0.02, 0.01],
        [0.22, 0.28, 0.2, 0.04, -0.05, 0.0, 0.01, 0.01, -0.01],
        [0.21, 0.20, 0.28, 0.05, -0.06, 0.0, 0.04, 0.03, 0.01],
        [0.09, 0.04, 0.05, 0.26, 0.14, 0.05, 0.05, 0.05, 0.04],
        [-0.02, -0.05, -0.06, 0.14, 0.20, 0.07, 0.05, 0.05, 0.05],
        [0.01, 0.0, 0.0, 0.05, 0.07, 0.12, 0.08, 0.07, 0.06],
        [0.03, 0.01, 0.04, 0.05, 0.05, 0.08, 0.12, 0.1, 0.08],
        [0.02, 0.01, 0.03, 0.05, 0.05, 0.07, 0.1, 0.1, 0.09],
        [0.01, -0.01, 0.01, 0.04, 0.05, 0.06, 0.08, 0.09, 0.09]
    ])

    # Find the interval in which input period is located
    for i in range(len(Tlist) - 1):
        if T1 - Tlist[i] >= 0 and T1 - Tlist[i + 1] < 0:
            index1 = i
        if T2 - Tlist[i] >= 0 and T2 - Tlist[i + 1] < 0:
            index2 = i

    # Linearly interpolate the corresponding value of each coregionalization matrix coefficient
    def interpolate(B, index1, index2, T1, T2, Tlist):
        coeff1 = B[index1, index2] + (B[index1 + 1, index2] - B[index1, index2]) / (
            Tlist[index1 + 1] - Tlist[index1]) * (T1 - Tlist[index1])
        coeff2 = B[index1, index2 + 1] + (B[index1 + 1, index2 + 1] - B[index1, index2 + 1]) / (
            Tlist[index1 + 1] - Tlist[index1]) * (T1 - Tlist[index1])
        coeff0 = coeff1 + (coeff2 - coeff1) / \
            (Tlist[index2 + 1] - Tlist[index2]) * (T2 - Tlist[index2])
        return coeff0

    B1coeff0 = interpolate(B1, index1, index2, T1, T2, Tlist)
    B2coeff0 = interpolate(B2, index1, index2, T1, T2, Tlist)
    B3coeff0 = interpolate(B3, index1, index2, T1, T2, Tlist)

    # Compute the correlation coefficient (Equation 42)
    rho = B1coeff0 * np.exp(-3 * h / 20) + B2coeff0 * np.exp(-3 * h / 70)

    if h == 0:
        rho = B1coeff0 * np.exp(-3 * h / 20) + B2coeff0 * \
            np.exp(-3 * h / 70) + B3coeff0

    return rho

