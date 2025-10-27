"""
Created by Vitor Monteiro, 23/09/2025

Compute cross-spatial correlation coefficients for Sa(T) for
NGA-West2 for periods between [0.01-5.0]s.

For more details please see: 
Markhvida, M., Ceferino, L., & Baker, J. W. (2018). 
Modeling spatially correlated spectral accelerations at multiple periods using principal
component analysis and geostatistics. Earthquake Engineering & Structural Dynamics, 47(5), 1107–1123.
https://doi.org/10.1002/eqe.3007
"""


import numpy as np
from scipy.interpolate import interp1d

T_list = np.array([0.01, 0.02, 0.03, 0.05, 0.075,
                   0.1, 0.15, 0.2, 0.25, 0.3, 0.4,
                   0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0])

pcs = np.array([
    [0.27, -0.14, 0.07, -0.11, -0.09],
    [0.27, -0.14, 0.08, -0.12, -0.10],
    [0.27, -0.15, 0.10, -0.14, -0.13],
    [0.25, -0.18, 0.18, -0.22, -0.18],
    [0.24, -0.22, 0.24, -0.23, -0.13],
    [0.23, -0.23, 0.23, -0.16, 0.04],
    [0.24, -0.21, 0.13, 0.08, 0.33],
    [0.25, -0.17, -0.01, 0.28, 0.40],
    [0.25, -0.12, -0.15, 0.37, 0.25],
    [0.25, -0.07, -0.24, 0.36, 0.04],
    [0.25, 0.01, -0.33, 0.23, -0.26],
    [0.25, 0.08, -0.36, 0.06, -0.34],
    [0.23, 0.19, -0.34, -0.22, -0.17],
    [0.21, 0.26, -0.24, -0.33, 0.08],
    [0.19, 0.33, -0.09, -0.27, 0.36],
    [0.18, 0.36, 0.06, -0.16, 0.35],
    [0.17, 0.36, 0.26, 0.07, 0.06],
    [0.16, 0.35, 0.35, 0.24, -0.16],
    [0.15, 0.33, 0.37, 0.33, -0.28]
])

nested_para = np.array([
    [2.50, 4.52, 15, 6.78, 250],
    [0.50, 1.40, 10, 2.60, 160],
    [0.15, 0.42, 15, 0.63, 160],
    [0.15, 0.23, 10, 0.23, 120],
    [0.31, 0.0001, 0.0001, 0.0001, 0.0001],
])

# Interpolators for PCs
pcs_interp = [interp1d(T_list, pcs[:, i], kind='linear',
                       fill_value="extrapolate") for i in range(5)]


def get_pc(T):
    return np.array([f(T) for f in pcs_interp])


def CrossSpatialCorrMCB18(T1, T2, h):
    # Interpolated PCs
    pc1 = get_pc(T1)
    pc2 = get_pc(T2)

    # Nugget
    Inugget = 1 if h == 0 else 0

    # Cross-covariance
    C_h = 0.0
    Cii_0 = 0.0
    Cjj_0 = 0.0
    for i in range(5):
        c0i, c1i, a1i, c2i, a2i = nested_para[i]
        Cij_h = (c0i * (Inugget) + c1i * (np.exp(-3 * h / a1i)) +
                 c2i * (np.exp(-3 * h / a2i)))
        # Total sill component
        silli = c0i + c1i + c2i

        C_h += pc1[i] * pc2[i] * (Cij_h) / 0.95

        Cii_0 += pc1[i] * pc1[i] * silli / 0.95
        Cjj_0 += pc2[i] * pc2[i] * silli / 0.95

    C_h_other = (C_h / np.sqrt(Cii_0 * Cjj_0))

    return C_h_other

