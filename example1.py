import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

from models.monteiroEtAl26 import CrossSpatialCorrMAO26
from models.markhvidaEtAl18 import CrossSpatialCorrMCB18
from models.lothbaker13 import CrossSpatialCorrLB13
from models.duning21 import CrossSpatialCorrDN21
from models.jayarambaker09 import SpatialCorrJB09

# --- Setup ---
my_path = Path(__file__).parent
chosen_period = 0.1
periods = np.linspace(0.2 * chosen_period, 2.0 * chosen_period, 10)
h_model = np.linspace(0, 150, 151)

# --- Load indirect correlation data ---
file_path = my_path / f'outputs_within/WithinCorr_Saavg2_ind/WithinCorrSaavg2({chosen_period:.2f})ind.csv'
cols = ['Bin', 'Correlation_loth', 'Correlation_markhvida', 'Correlation_DuNing', 'Correlation_vitor', 'Correlation_markov']
df = pd.read_csv(file_path, usecols=cols)

# --- Define model functions and colors ---
models = {
    'LB13':   (CrossSpatialCorrLB13,    'lightcoral'),
    'MCB18':  (CrossSpatialCorrMCB18,   'lightgreen'),
    'DN21':   (CrossSpatialCorrDN21,    'lightblue'),
    'MAO26':  (CrossSpatialCorrMAO26,   'lightsalmon'),
    'Markov': (lambda p1, p2, h: SpatialCorrJB09(p1, h, 1), 'lightgray')
}

# --- Compute direct MAO26 correlation for chosen period ---
im = f"Saavg2({chosen_period})"
direct_mao = [CrossSpatialCorrMAO26(im, im, h, cluster=0) for h in h_model]

# --- Compute indirect (Sa(T)) curves for each model ---
model_curves = {
    name: [[
        (func(f"Sa({t})", f"Sa({t})", h, cluster=0) if name == "MAO26"
         else func(t, t, h) if name != "DN21"
         else func(f"SA({t})", f"SA({t})", h))
        for h in h_model]
        for t in periods]
    for name, (func, _) in models.items()
}

# --- Plot ---
plt.figure(figsize=(8, 6))

# Plot background Sa(T) curves
for name, (_, color) in models.items():
    for curve in model_curves[name]:
        plt.plot(h_model, curve, color=color, lw=1, alpha=0.6)

plt.plot([], [], color='lightgray', linewidth=1, label='Sa(T) curves', alpha=0.6)
        

# Plot indirect correlations
plt.plot(df['Bin'], df['Correlation_loth'], label='Indirect - LB13', color='r', lw=2)
plt.plot(df['Bin'], df['Correlation_markhvida'], label='Indirect - MCB18', color='g', lw=2)
plt.plot(df['Bin'], df['Correlation_DuNing'], label='Indirect - DN21', color='b', lw=2)
plt.plot(df['Bin'], df['Correlation_vitor'], label='Indirect - MAO25', color='salmon', lw=2)
plt.plot(df['Bin'], df['Correlation_markov'], label='Indirect - Markov-hypothesis', color='gray', lw=2)

# Plot direct curve
plt.plot(h_model, direct_mao, label='Direct - MAO26', color='k', lw=2)

# --- Aesthetics ---
plt.xlabel('Distance, h [km]', fontsize=18)
plt.ylabel(fr'$\rho_{{\ln Sa_{{avg2}}({chosen_period})_n,\ln Sa_{{avg2}}({chosen_period})_m}}$', fontsize=26)
plt.xlim(0, 150); plt.ylim(0, 1)
plt.xticks(fontsize=14); plt.yticks(fontsize=14)
plt.grid(True)
plt.legend(fontsize=14)
plt.tight_layout()


# Save or show
figures_folder = my_path / "Figures"
figures_folder.mkdir(parents=True, exist_ok=True)
plt.savefig(figures_folder / f"comp_avgsa2_sa_dir_indir_{chosen_period}.pdf",
            dpi=300, bbox_inches='tight')
plt.show()
