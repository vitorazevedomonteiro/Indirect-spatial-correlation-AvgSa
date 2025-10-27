from pathlib import Path
import pandas as pd
import numpy as np
import itertools

# Import spatial models
from models.lothbaker13 import CrossSpatialCorrLB13
from models.markhvidaEtAl18 import CrossSpatialCorrMCB18
from models.duning21 import CrossSpatialCorrDN21
from models.monteiroEtAl26 import MAO2026
from models.jayarambaker09 import SpatialCorrJB09

# Import non-spatial models
from models.bakerjayaram08 import corrBJ08



import sys
import os


# Define paths
my_path = Path(__file__).parent

# Periods to analyze
avgsa_periods = [0.10, 0.15, 0.20, 0.30, 0.40,
                 0.50, 0.60, 0.75, 0.80, 0.90, 1.00, 1.20, 1.50, 2.0, 2.5]

# Function to process and create the new CSV with the desired columns

def process_period(period):
    # Load the standard deviation data for the current period
    stdev_sa_path = my_path / f"predicted_files/predicted_Sa_avg2_sa({period:.2f}).csv"
    stdev_sa_df = pd.read_csv(stdev_sa_path)

    # Get the data for the first RSN only
    first_rsn = stdev_sa_df['RSN'].iloc[0]
    first_rsn_data = stdev_sa_df[stdev_sa_df['RSN'] == first_rsn]

    # Map Period to Stdev3
    period_to_stdev = dict(
        zip(first_rsn_data['Period'], first_rsn_data['Stdev3']))
    periods = list(period_to_stdev.keys())

    # Cartesian product (all possible pairs including reverse and self-pairs)
    all_combos = itertools.product(periods, periods)

    rows = []
    for p1, p2 in all_combos:
        rows.append({
            'Period1': p1,
            'Period2': p2,
            'stdev3_period1': period_to_stdev[p1],
            'stdev3_period2': period_to_stdev[p2]
        })

    # Save result
    output_folder = my_path / "stdev_combinations"
    output_folder.mkdir(parents=True, exist_ok=True)

    df_result = pd.DataFrame(rows)
    df_result.to_csv(
        output_folder / f"stdev_combinations_{period:.2f}.csv", index=False)


# Process each period in avgsa_periods
for period in avgsa_periods:
    process_period(period)


# After processing stdev combinations, now compute correlations
for period in avgsa_periods:
    # Load the stdev combination file
    stdev_file = my_path / \
        f"stdev_combinations/stdev_combinations_{period:.2f}.csv"
    stdev_df = pd.read_csv(stdev_file)

    bin_values = np.linspace(0, 150, 151)
    correlation_rows = []

    for _, row in stdev_df.iterrows():
        p1 = row['Period1']
        p2 = row['Period2']
        s1 = row['stdev3_period1']
        s2 = row['stdev3_period2']

        for bin_val in bin_values:
            rho_loth = CrossSpatialCorrLB13(p1, p2, bin_val)
            rho_markhvida = CrossSpatialCorrMCB18(p1, p2, bin_val)[3]
            rho_DuNing = CrossSpatialCorrDN21(f"SA({p1})", f"SA({p2})", bin_val)[3]
            rho_vitor = MAO2026(
                f"SA({p1})", f"SA({p2})", bin_val, "NGA_ESM_database2")
            rho_markov = corrBJ08(p1, p2) * SpatialCorrJB09(np.max([p1, p2]), bin_val, 1)
            correlation_rows.append({
                'Bin': bin_val,
                'Period1': p1,
                'Period2': p2,
                'stdev3_period1': s1,
                'stdev3_period2': s2,
                'Correlation_loth': rho_loth,
                'Correlation_markhvida': rho_markhvida,
                'Correlation_DuNing': rho_DuNing,
                'Correlation_vitor': rho_vitor,
                'Correlation_markov': rho_markov,
                'numerator_loth': rho_loth * s1 * s2,
                'numerator_markhvida': rho_markhvida * s1 * s2,
                'numerator_DuNing': rho_DuNing * s1 * s2,
                'numerator_vitor': rho_vitor * s1 * s2,
                'numerator_markov': rho_markov * s1 * s2,
            })

    # Save the correlation file
    correlation_df = pd.DataFrame(correlation_rows)
    correlation_df.to_csv(
        my_path / f"stdev_combinations/correlation_sa_{period:.2f}.csv", index=False)



# After saving correlation files, compute the summed numerator per bin and save
for period in avgsa_periods:
    correlation_path = my_path / \
        f"stdev_combinations/correlation_sa_{period:.2f}.csv"
    correlation_df = pd.read_csv(correlation_path)

    # Group by Bin, sum numerators separately
    numerator_loth = correlation_df.groupby("Bin", as_index=False)[
        "numerator_loth"].sum()
    numerator_markhvida = correlation_df.groupby("Bin", as_index=False)[
        "numerator_markhvida"].sum()
    numerator_DuNing = correlation_df.groupby("Bin", as_index=False)[
        "numerator_DuNing"].sum()
    numerator_vitor = correlation_df.groupby("Bin", as_index=False)[
        "numerator_vitor"].sum()
    numerator_markov = correlation_df.groupby("Bin", as_index=False)[
        "numerator_markov"].sum()

    # Merge the four results on Bin, pairwise
    numerator_df = pd.merge(numerator_loth, numerator_markhvida, on="Bin")
    numerator_df = pd.merge(numerator_df, numerator_DuNing, on="Bin")
    numerator_df = pd.merge(numerator_df, numerator_vitor, on="Bin")
    numerator_df = pd.merge(numerator_df, numerator_markov, on="Bin")

    # Divide each by 100
    numerator_df["numerator_loth"] /= 100
    numerator_df["numerator_markhvida"] /= 100
    numerator_df["numerator_DuNing"] /= 100
    numerator_df["numerator_vitor"] /= 100
    numerator_df["numerator_markov"] /= 100

    # Save result
    numerator_df.to_csv(
        my_path / f"stdev_combinations/numerator_{period:.2f}.csv", index=False
    )


# After processing stdev combinations, now compute zero-distance correlations and denominators
for period in avgsa_periods:
    # Load the stdev combination file
    stdev_file = my_path / \
        f"stdev_combinations/stdev_combinations_{period:.2f}.csv"
    stdev_df = pd.read_csv(stdev_file)

    correlation_rows = []

    for _, row in stdev_df.iterrows():
        p1 = row['Period1']
        p2 = row['Period2']
        s1 = row['stdev3_period1']
        s2 = row['stdev3_period2']

        rho_loth = CrossSpatialCorrLB13(p1, p2, 0)
        rho_markhvida = CrossSpatialCorrMCB18(p1, p2, 0)[3]
        rho_DuNing = CrossSpatialCorrDN21(f"SA({p1})", f"SA({p2})", 0)[3]
        rho_vitor = MAO2026(f"SA({p1})", f"SA({p2})", 0, "NGA_ESM_database2")
        rho_markov = corrBJ08(p1,p2)
        correlation_rows.append({
            'Period1': p1,
            'Period2': p2,
            'stdev3_period1': s1,
            'stdev3_period2': s2,
            'Correlation_loth': rho_loth,
            'Correlation_markhvida': rho_markhvida,
            'Correlation_DuNing': rho_DuNing,
            'Correlation_vitor': rho_vitor,
            'Correlation_markov': rho_markov,
            'denominator_loth': rho_loth * s1 * s2,
            'denominator_markhvida': rho_markhvida * s1 * s2,
            'denominator_DuNing': rho_DuNing * s1 * s2,
            'denominator_vitor': rho_vitor * s1 * s2,
            'denominator_markov': rho_markov * s1 * s2
        })

    # Create correlation DataFrame
    correlation_df = pd.DataFrame(correlation_rows)

    # Sum all denominators and divide by 100
    total_denominator_loth = correlation_df['denominator_loth'].sum() / 100
    total_denominator_markhvida = correlation_df['denominator_markhvida'].sum(
    ) / 100
    total_denominator_DuNing = correlation_df['denominator_DuNing'].sum(
    ) / 100
    total_denominator_vitor = correlation_df['denominator_vitor'].sum(
    ) / 100
    total_denominator_markov = correlation_df['denominator_markov'].sum(
    ) / 100

    # Create a one-row DataFrame for output
    denominator_df = pd.DataFrame([{
        'denominator_loth': total_denominator_loth,
        'denominator_markhvida': total_denominator_markhvida,
        'denominator_DuNing': total_denominator_DuNing,
        'denominator_vitor': total_denominator_vitor,
        'denominator_markov': total_denominator_markov
    }])

    # Save the result
    denominator_df.to_csv(
        my_path / f"stdev_combinations/denominator_{period:.2f}.csv", index=False
    )
