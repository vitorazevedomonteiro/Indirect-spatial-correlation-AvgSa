from pathlib import Path
import pandas as pd
import numpy as np
import itertools
import sys
import os

# Import spatial models
from models.lothbaker13 import CrossSpatialCorrLB13
from models.markhvidaEtAl18 import CrossSpatialCorrMCB18
from models.duning21 import CrossSpatialCorrDN21
from models.monteiroEtAl26 import CrossSpatialCorrMAO26
from openquake.hazardlib.cross_correlation import GodaAtkinson2009
from openquake.hazardlib.imt import SA
ga09_model = GodaAtkinson2009()

# Import non-spatial models
from models.bakerjayaram08 import corrBJ08


def process_stdev_combinations_total(predicted_dir: Path, output_dir: Path, avgsa_periods):
    """
    Generates stdev combinations (Period1, Period2, stdev3_period1, stdev3_period2)
    from GMM-predicted files.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    for period in avgsa_periods:
        file_path = predicted_dir / f"predicted_Saavg2_sa({period:.2f}).csv"
        df = pd.read_csv(file_path)

        # Use first RSN only (since all RSNs share same Period-Stdev structure)
        first_rsn = df['RSN'].iloc[0]
        first_rsn_data = df[df['RSN'] == first_rsn]

        period_to_stdev1 = dict(zip(first_rsn_data['Period'], first_rsn_data['Stdev1']))
        period_to_stdev2 = dict(zip(first_rsn_data['Period'], first_rsn_data['Stdev2']))
        period_to_stdev3 = dict(zip(first_rsn_data['Period'], first_rsn_data['Stdev3']))
        periods = list(period_to_stdev3.keys())

        combos = itertools.product(periods, periods)
        rows = [
            {'Period1': p1, 'Period2': p2,
            'stdev1_period1': period_to_stdev1[p1],
            'stdev1_period2': period_to_stdev1[p2],
            'stdev2_period1': period_to_stdev2[p1],
            'stdev2_period2': period_to_stdev2[p2],
            'stdev3_period1': period_to_stdev3[p1],
            'stdev3_period2': period_to_stdev3[p2]}
            for p1, p2 in combos
        ]

        output_file = output_dir / f"stdev_combinations_{period:.2f}.csv"
        pd.DataFrame(rows).to_csv(output_file, index=False)
        print(f"Created: {output_file.name}")


def compute_correlations_total(stdev_dir: Path, avgsa_periods):
    """
    Computes correlation and numerator terms across distance bins
    for all spatial/non-spatial models.
    """
    for period in avgsa_periods:
        stdev_file = stdev_dir / f"stdev_combinations_{period:.2f}.csv"
        stdev_df = pd.read_csv(stdev_file)

        bin_values = np.linspace(0, 150, 151)
        correlation_rows = []

        for _, row in stdev_df.iterrows():
            p1, p2 = row['Period1'], row['Period2']
            std2_1, std2_2 = row['stdev2_period1'], row['stdev2_period2']
            std3_1, std3_2 = row['stdev3_period1'], row['stdev3_period2']

            for bin_val in bin_values:
                rho_loth = CrossSpatialCorrLB13(p1, p2, bin_val)
                rho_markhvida = CrossSpatialCorrMCB18(p1, p2, bin_val)
                rho_DuNing = CrossSpatialCorrDN21(f"SA({p1})", f"SA({p2})", bin_val)
                rho_vitor = CrossSpatialCorrMAO26(f"Sa({p1})", f"Sa({p2})", bin_val, cluster=0)
                imt1 = SA(p1)
                imt2 = SA(p2)
                rho_between = ga09_model.get_correlation(imt1, imt2)

                correlation_rows.append({
                    'Bin': bin_val,
                    'Period1': p1,
                    'Period2': p2,
                    'stdev2_period1': std2_1,
                    'stdev2_period2': std2_2,
                    'stdev3_period1': std3_1,
                    'stdev3_period2': std3_2,
                    'Correlation_loth': rho_loth,
                    'Correlation_markhvida': rho_markhvida,
                    'Correlation_DuNing': rho_DuNing,
                    'Correlation_vitor': rho_vitor,
                    'Correlation_between': rho_between,
                    'numerator_between': rho_between * std2_1 * std2_2,
                    'numerator_between': rho_between * std2_1 * std2_2,
                    'numerator_loth': (rho_between * std2_1 * std2_2 + rho_loth * std3_1 * std3_2),
                    'numerator_markhvida': (rho_between * std2_1 * std2_2 + rho_markhvida * std3_1 * std3_2),
                    'numerator_DuNing': (rho_between * std2_1 * std2_2 + rho_DuNing * std3_1 * std3_2),
                    'numerator_vitor': (rho_between * std2_1 * std2_2 + rho_vitor * std3_1 * std3_2),
                })

        output_file = stdev_dir / f"correlation_sa_{period:.2f}.csv"
        pd.DataFrame(correlation_rows).to_csv(output_file, index=False)
        print(f"Correlation file saved: {output_file.name}")


def compute_numerators_total(stdev_dir: Path, avgsa_periods):
    """
    Aggregates numerators by Bin for each model and saves as separate CSVs.
    """
    for period in avgsa_periods:
        correlation_path = stdev_dir / f"correlation_sa_{period:.2f}.csv"
        correlation_df = pd.read_csv(correlation_path)

        grouped = {
            name: correlation_df.groupby("Bin", as_index=False)[name].sum()
            for name in [
                "numerator_loth", "numerator_markhvida",
                "numerator_DuNing", "numerator_vitor"
            ]
        }

        numerator_df = grouped["numerator_loth"]
        for key, df in grouped.items():
            if key != "numerator_loth":
                numerator_df = numerator_df.merge(df, on="Bin")

        # Scale down (as in your original script)
        for col in numerator_df.columns:
            if col.startswith("numerator_"):
                numerator_df[col] /= 100

        output_file = stdev_dir / f"numerator_{period:.2f}.csv"
        numerator_df.to_csv(output_file, index=False)
        print(f"Numerator saved: {output_file.name}")


def compute_denominators_total(stdev_dir: Path, avgsa_periods):
    """
    Computes denominators (zero-distance correlations) for all models.
    """
    for period in avgsa_periods:
        stdev_file = stdev_dir / f"stdev_combinations_{period:.2f}.csv"
        stdev_df = pd.read_csv(stdev_file)

        rows = []
        for _, row in stdev_df.iterrows():
            p1, p2 = row['Period1'], row['Period2']
            std2_1, std2_2 = row['stdev2_period1'], row['stdev2_period2']
            std3_1, std3_2 = row['stdev3_period1'], row['stdev3_period2']

            rho_loth = CrossSpatialCorrLB13(p1, p2, 0)
            rho_markhvida = CrossSpatialCorrMCB18(p1, p2, 0)
            rho_DuNing = CrossSpatialCorrDN21(f"SA({p1})", f"SA({p2})", 0)
            rho_vitor = CrossSpatialCorrMAO26(f"Sa({p1})", f"Sa({p2})", 0, cluster=0)
            imt1 = SA(p1)
            imt2 = SA(p2)
            rho_between = ga09_model.get_correlation(imt1, imt2)
        
            rows.append({
                'Correlation_between': rho_between,
                'denominator_loth': (rho_between * std2_1 * std2_2 + rho_loth * std3_1 * std3_2),
                'denominator_markhvida': (rho_between * std2_1 * std2_2 + rho_markhvida * std3_1 * std3_2),
                'denominator_DuNing': (rho_between * std2_1 * std2_2 + rho_DuNing * std3_1 * std3_2),
                'denominator_vitor': (rho_between * std2_1 * std2_2 + rho_vitor * std3_1 * std3_2),
            })

        df = pd.DataFrame(rows)

        denominator_df = pd.DataFrame([{
            f"{key}": df[key].sum() / 100
            for key in df.columns
        }])

        output_file = stdev_dir / f"denominator_{period:.2f}.csv"
        denominator_df.to_csv(output_file, index=False)
        print(f"Denominator saved: {output_file.name}")