from pathlib import Path
from utils.gmm_calculator import compute_gmm_predictions
from utils.stdev_and_corr import (
    process_stdev_combinations, compute_correlations,
    compute_numerators, compute_denominators
)
from utils.correlation_final import compute_final_correlations

def main():
    base_path = Path(__file__).parent
    data_path = base_path / "data" / "Database.csv"
    predicted_dir = base_path / "outputs" / "predicted_files"
    stdev_dir = base_path / "outputs" / "stdev_comb"
    corr_output_dir = base_path / "outputs" / "Corr_AvgSa2_Ind"

    
    avgsa_periods = [0.1, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.75, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0, 2.5]

    print("\n=== STEP 1: GMM PREDICTIONS ===")
    compute_gmm_predictions(data_path, predicted_dir, avgsa_periods)

    print("\n=== STEP 2: STDEV COMBINATIONS & CORRELATIONS ===")
    process_stdev_combinations(predicted_dir, stdev_dir, avgsa_periods)
    compute_correlations(stdev_dir, avgsa_periods)
    compute_numerators(stdev_dir, avgsa_periods)
    compute_denominators(stdev_dir, avgsa_periods)

    print("\n=== STEP 3: FINAL CORRELATIONS ===")
    compute_final_correlations(stdev_dir, corr_output_dir, avgsa_periods)

    print("\nAll processing complete!")

if __name__ == "__main__":
    main()
