"""
Created by Vitor Monteiro, 23/09/2025

Compute spatial correlation coefficients for Sa(T) and PGA for
NGA-W1 database for periods between [0.01-10]s.

For more details please see: 
Jayaram, N., & Baker, J. W. (2009). Correlation model for
spatially distributed ground‐motion intensities. 
Earthquake Engineering & Structural Dynamics, 38(15), 1687–1708.
https://doi.org/10.1002/eqe.922
"""


import numpy as np

def SpatialCorrJB09(T, h, vs30=1):
    # Use T=0 for PGA
    if not (0 <= T <= 10.0):
        raise ValueError(f"T = {T} is outside the valid range [0, 10.0].")
        
    if vs30 == 1:  # vs30 clustering 
        if T < 1:
            b = 8.5 + 17.2 * T
        else:  # If vs30 clustering
            b = 22.0 + 3.7 * T
    elif vs30 == 2:
        if T < 1:
            b = 40.7 - 15.0 * T
        else:
            b = 22.0 + 3.7 * T
    return np.exp(-(3 * h) / b)


