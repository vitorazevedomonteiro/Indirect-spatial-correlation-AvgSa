"""
Created by Vitor Monteiro, 23/09/2025

Compute cross-spatial correlation coefficients for Sa(T), PGA, PGV, and Duration for
NGA-West2 for periods between [0.01-10.0]s.

For more details please see: 
Du, W., & Ning, C.-L. (2021). Modeling spatial cross-correlation of
multiple ground motion intensity measures (SAs, PGA, PGV, Ia, CAV, and
significant durations) based on principal component and geostatistical analyses.
Earthquake Spectra, 37(1), 486–504. https://doi.org/10.1177/8755293020952442
"""


import numpy as np
from scipy.interpolate import interp1d

# Use numeric placeholders for PGA and PGV: -1.0 for PGA, -2.0 for PGV
T_list = np.array(['SA(0.01)', 'SA(0.05)', 'SA(0.075)',
                   'SA(0.1)', 'SA(0.2)', 'SA(0.3)', 'SA(0.4)',
                   'SA(0.5)', 'SA(0.75)', 'SA(1.0)', 'SA(1.5)',
                   'SA(2.0)', 'SA(3.0)', 'SA(4.0)', 'SA(5.0)',
                   'SA(7.5)', 'SA(10)', 'PGA', 'PGV', 'DS575', 'DS595'])

pcs = np.array([
    [0.28, -0.15, 0.07, 0.04, -0.05, -0.10, -0.08],
    [0.24, -0.19, 0.23, 0.04, -0.18, -0.10, 0.06],
    [0.21, -0.21, 0.29, 0.08, -0.17, 0.06, 0.34],
    [0.20, -0.22, 0.27, 0.08, -0.09, 0.15, 0.45],
    [0.22, -0.20, 0.09, 0.04, 0.23, 0.30, 0.11],
    [0.24, -0.14, -0.11, 0.02, 0.39, 0.28, -0.11],
    [0.24, -0.09, -0.19, -0.05, 0.37, 0.20, -0.13],
    [0.24, -0.04, -0.25, -0.11, 0.31, 0.02, -0.01],
    [0.23, 0.05, -0.30, -0.20, 0.08, -0.19, 0.22],
    [0.22, 0.13, -0.27, -0.27, -0.10, -0.17, 0.32],
    [0.19, 0.22, -0.20, -0.26, -0.23, 0.02, 0.19],
    [0.19, 0.27, -0.10, -0.17, -0.27, 0.20, 0.06],
    [0.17, 0.30, 0.06, 0.00, -0.25, 0.38, -0.20],
    [0.16, 0.32, 0.12, 0.11, -0.11, 0.32, -0.14],
    [0.15, 0.32, 0.15, 0.24, 0.06, 0.17, -0.08],
    [0.12, 0.31, 0.08, 0.37, 0.23, -0.25, 0.16],
    [0.12, 0.30, 0.01, 0.41, 0.20, -0.30, 0.18],
    [0.28, -0.15, 0.07, 0.04, -0.05, -0.10, -0.08],
    [0.26, 0.13, -0.08, 0.09, -0.10, -0.29, -0.23],
    [-0.09, 0.21, 0.41, -0.40, 0.29, -0.09, 0.10],
    [-0.07, 0.27, 0.36, -0.38, 0.27, 0.02, 0.09]
])

nested_para = np.array([
    [1.03, 0.88, 15, 10.11, 200],
    [0.36, 1.76, 25, 2.61, 150],
    [0.13, 0.37, 25, 1.75, 200],
    [0.09, 0.26, 20, 1.11, 150],
    [0.10, 0.32, 15, 0.45, 150],
    [0.11, 0.13, 25, 0.35, 250],
    [0.06, 0.16, 25, 0.25, 250]
])

# Create a dictionary to quickly fetch PC values using label
pcs_dict = {label: pcs[i, :] for i, label in enumerate(T_list)}


def get_pc(label):
    """Fetch or interpolate principal component values by IM label."""
    if label in pcs_dict:
        return pcs_dict[label]

    # Check if label looks like 'SA(some_number)'
    if label.startswith("SA(") and label.endswith(")"):
        try:
            period = float(label[3:-1])
        except ValueError:
            raise ValueError(f"Could not parse period from label '{label}'.")

        # Make sure period is within interpolation range
        if not (0.01 <= period <= 10):
            raise ValueError(
                f"Period {period} is outside the interpolation range [0.01, 10].")

        # Extract known SA periods and PCs
        sa_periods = []
        sa_pcs = []

        for i, name in enumerate(T_list):
            if name.startswith("SA("):
                try:
                    p = float(name[3:-1])
                    sa_periods.append(p)
                    sa_pcs.append(pcs[i, :])
                except ValueError:
                    continue

        sa_periods = np.array(sa_periods)
        sa_pcs = np.array(sa_pcs)

        # Interpolate each PC component in log-space
        log_periods = np.log(sa_periods)
        log_target = np.log(period)

        interpolated_pc = np.zeros(pcs.shape[1])
        for i in range(pcs.shape[1]):
            f_interp = interp1d(
                log_periods, sa_pcs[:, i], kind='linear', fill_value='extrapolate')
            interpolated_pc[i] = f_interp(log_target)
        return interpolated_pc

    raise ValueError(
        f"IM label '{label}' not found and cannot be interpolated.")




def CrossSpatialCorrDN21(IM1, IM2, h):
    # Interpolated PCs
    pc1 = get_pc(IM1)
    pc2 = get_pc(IM2)

    Inugget = 1 if h == 0 else 0
    
    C_h = 0.0
    Cii_0 = 0.0
    Cjj_0 = 0.0
    for i in range(7):
        c0i, c1i, a1i, c2i, a2i = nested_para[i]
        gammai_h = (c0i * (1 - Inugget) + c1i * (1 - np.exp(-3 * h / a1i)) +
                    c2i * (1 - np.exp(-3 * h / a2i)))
        Cij_h = (c0i * (Inugget) + c1i * (np.exp(-3 * h / a1i)) +
                 c2i * (np.exp(-3 * h / a2i)))
        silli = c0i + c1i + c2i


        C_h += pc1[i] * pc2[i] * (Cij_h) / 0.90
        Cii_0 += pc1[i] * pc1[i] * silli / 0.9
        Cjj_0 += pc2[i] * pc2[i] * silli / 0.9

    C_h_other = (C_h / np.sqrt(Cii_0 * Cjj_0))

    return C_h_other

