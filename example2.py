import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
from pathlib import Path
from models.monteiroEtAl26 import CrossSpatialCorrMAO26
from models.markhvidaEtAl18 import CrossSpatialCorrMCB18
from models.lothbaker13 import CrossSpatialCorrLB13
from models.duning21 import CrossSpatialCorrDN21
from matplotlib.lines import Line2D

# --- Setup ---
my_path = Path(__file__).parent
chosen_period = 1.0
h_model = np.linspace(0, 150, 151)

# --- Load indirect correlation data ---
file_path = my_path / f'outputs_total/TotalCorr_Saavg2_ind/TotalCorrSaavg2({chosen_period:.2f})ind.csv'
cols = ['Bin', 'Correlation_loth', 'Correlation_markhvida', 'Correlation_DuNing', 'Correlation_vitor']
df = pd.read_csv(file_path, usecols=cols)

# --- Load indirect correlation data ---
file_path_w = my_path / f'outputs_within/WithinCorr_Saavg2_ind/WithinCorrSaavg2({chosen_period:.2f})ind.csv'
cols_w = ['Bin', 'Correlation_loth', 'Correlation_markhvida', 'Correlation_DuNing', 'Correlation_vitor']
df_w = pd.read_csv(file_path_w, usecols=cols)


plt.figure(figsize=(8, 6))

# This file contains spatial correlations using the indirect approach with total-event residuals 
plt.plot(df['Bin'], df['Correlation_loth'], color='r',linestyle='--', lw=2)
plt.plot(df['Bin'], df['Correlation_markhvida'], color='g', linestyle='--', lw=2)
plt.plot(df['Bin'], df['Correlation_DuNing'], color='b', linestyle='--', lw=2)
plt.plot(df['Bin'], df['Correlation_vitor'], color='orange', linestyle='--', lw=2)

# This file contains spatial correlations using the indirect approach with within-event residuals
plt.plot(df_w['Bin'], df_w['Correlation_loth'], label='Indirect - LB13', color='r', lw=2)
plt.plot(df_w['Bin'], df_w['Correlation_markhvida'], label='Indirect - MCB18', color='g', lw=2)
plt.plot(df_w['Bin'], df_w['Correlation_DuNing'], label='Indirect - DN21', color='b', lw=2)
plt.plot(df_w['Bin'], df_w['Correlation_vitor'], label='Indirect - MAO26', color='orange', lw=2)


plt.xlabel('h [km]', fontsize=18)
plt.xlim(0, 150)
plt.ylim(0, 1)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

plt.ylabel(fr'$\rho_{{\ln Sa_{{avg2}}({chosen_period})_n,\ln Sa_{{avg2}}({chosen_period})_m}}$', fontsize=26)
plt.title('')
plt.grid(True)


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
plt.savefig(figures_folder / f"comp_avgsa2_indir_within_total_{chosen_period:.2f}.png",
            dpi=300, bbox_inches='tight')
plt.show()
