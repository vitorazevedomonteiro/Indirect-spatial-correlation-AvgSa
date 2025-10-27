from pathlib import Path
from utils.gmm_calculator import compute_gmm_predictions
from utils.stdev_and_corr import (
    process_stdev_combinations, compute_correlations,
    compute_numerators, compute_denominators)
from utils.corr_final import compute_final_corr
from utils.stdev_and_corr_total import (
    process_stdev_combinations_total, compute_correlations_total,
    compute_numerators_total, compute_denominators_total)
from utils.corr_final_total import compute_final_corr_total



def main_within():
    im = 'Saavg2' # could be 'Saavg2' or 'Saavg3
    base_path = Path(__file__).parent
    data_path = base_path / "data" / "Database.csv"
    predicted_dir = base_path / "outputs_within" / "predicted_files"
    stdev_dir = base_path / "outputs_within" / "stdev_comb"
    corr_output_dir = base_path / "outputs_within" / f"WithinCorr_{im}_ind"


    avgsa_periods = [0.1, 1.0] 
    # add the periods that you want, but be aware of the GMM period limits!

    print("\n=== STEP 1: GMM PREDICTIONS ===")
    compute_gmm_predictions(im, data_path, predicted_dir, avgsa_periods)

    print("\n=== STEP 2: STDEV COMBINATIONS & CORRELATIONS ===")
    process_stdev_combinations(predicted_dir, stdev_dir, avgsa_periods)
    compute_correlations(stdev_dir, avgsa_periods)
    compute_numerators(stdev_dir, avgsa_periods)
    compute_denominators(stdev_dir, avgsa_periods)

    print("\n=== STEP 3: FINAL CORRELATIONS ===")
    compute_final_corr(stdev_dir, corr_output_dir, avgsa_periods)

    print("\nAll processing complete!")
    

def main_total():
    im = 'Saavg2' # could be 'Saavg2' or 'Saavg3
    base_path = Path(__file__).parent
    data_path = base_path / "data" / "Database.csv"
    predicted_dir = base_path / "outputs_total" / "predicted_files"
    stdev_dir = base_path / "outputs_total" / "stdev_comb"
    corr_output_dir = base_path / "outputs_total" / f"TotalCorr_{im}_ind"


    avgsa_periods = [1.0] 
    # add the periods that you want, but be aware of the GMM period limits!

    print("\n=== STEP 1: GMM PREDICTIONS ===")
    compute_gmm_predictions(im, data_path, predicted_dir, avgsa_periods)

    print("\n=== STEP 2: STDEV COMBINATIONS & CORRELATIONS ===")
    process_stdev_combinations_total(predicted_dir, stdev_dir, avgsa_periods)
    compute_correlations_total(stdev_dir, avgsa_periods)
    compute_numerators_total(stdev_dir, avgsa_periods)
    compute_denominators_total(stdev_dir, avgsa_periods)

    print("\n=== STEP 3: FINAL CORRELATIONS ===")
    compute_final_corr_total(stdev_dir, corr_output_dir, avgsa_periods)

    print("\nAll processing complete!")



if __name__ == "__main__":
    main_within()
