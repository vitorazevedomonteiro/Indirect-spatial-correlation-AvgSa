from pathlib import Path
import pandas as pd
import numpy as np


def compute_final_corr_total(stdev_dir: Path, output_dir: Path, avgsa_periods):
    """
    Combines numerator and denominator results to compute final correlations
    for all spatial/non-spatial models.

    Args:
        stdev_dir (Path): Directory containing numerator_*.csv and denominator_*.csv
        output_dir (Path): Output directory for final correlation CSVs
        avgsa_periods (list[float]): List of periods to process
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    for period in avgsa_periods:
        numerator_path = stdev_dir / f"numerator_{period:.2f}.csv"
        denominator_path = stdev_dir / f"denominator_{period:.2f}.csv"

        if not numerator_path.exists() or not denominator_path.exists():
            print(f"Missing files for period {period:.2f}, skipping.")
            continue

        numerator_df = pd.read_csv(numerator_path)
        denominator_df = pd.read_csv(denominator_path)

        # Extract single-row denominator values
        denominators = denominator_df.iloc[0].to_dict()

        # Build correlation DataFrame
        corr_df = pd.DataFrame()
        corr_df["Bin"] = numerator_df["Bin"]

        # Compute normalized correlation for each model
        for model in ["loth", "markhvida", "DuNing", "vitor"]:
            num_col = f"numerator_{model}"
            denom_val = denominators[f"denominator_{model}"]
            corr_df[f"Correlation_{model}"] = numerator_df[num_col] / (np.sqrt(denom_val) ** 2)

        # Save output
        output_file = output_dir / f"TotalCorrSaavg2({period:.2f})ind.csv"
        corr_df.to_csv(output_file, index=False)
        print(f"Final correlation saved: {output_file.name}")
