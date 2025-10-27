import pandas as pd
from pathlib import Path
import numpy as np

my_path = Path(__file__).parent

avgsa_periods = [0.1, 0.15, 0.20, 0.30, 0.40, 0.50,
                 0.60, 0.75, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0, 2.5]

output_folder = my_path / "Corr_AVGSA2_IND"
output_folder.mkdir(parents=True, exist_ok=True)

for i in avgsa_periods:
    files_path = my_path / "stdev_combinations"
    

    # Load numerator and denominator
    numerator_file = pd.read_csv(files_path / f"numerator_{i:.2f}.csv")
    denominator_file = pd.read_csv(files_path / f"denominator_{i:.2f}.csv")

    # Extract the denominator value (only one value in the file)
    denominator_value_loth = denominator_file['denominator_loth'].iloc[0]
    denominator_value_markhvida = denominator_file['denominator_markhvida'].iloc[0]
    denominator_value_DuNing = denominator_file['denominator_DuNing'].iloc[0]
    denominator_value_vitor = denominator_file['denominator_vitor'].iloc[0]
    denominator_value_markov = denominator_file['denominator_markov'].iloc[0]

    # Compute correlation
    corr_df = pd.DataFrame()
    corr_df["Bin"] = numerator_file["Bin"]
    corr_df["Correlation_loth"] = numerator_file['numerator_loth'] / \
        (np.sqrt(denominator_value_loth) * np.sqrt(denominator_value_loth))
    corr_df["Correlation_markhvida"] = numerator_file['numerator_markhvida'] / \
        (np.sqrt(denominator_value_markhvida)
         * np.sqrt(denominator_value_markhvida))
    corr_df["Correlation_DuNing"] = numerator_file['numerator_DuNing'] / \
        (np.sqrt(denominator_value_DuNing)
         * np.sqrt(denominator_value_DuNing))
    corr_df["Correlation_vitor"] = numerator_file['numerator_vitor'] / \
        (np.sqrt(denominator_value_vitor)
         * np.sqrt(denominator_value_vitor))
    corr_df["Correlation_markov"] = numerator_file['numerator_markov'] / \
        (np.sqrt(denominator_value_markov)
         * np.sqrt(denominator_value_markov))

    # Save to new file
    corr_df.to_csv(output_folder / f"Corr_AVGSA2_IND_{i:.2f}.csv", index=False)
