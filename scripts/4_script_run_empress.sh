#!/bin/bash
#############################################################
# Title: Automate Empress Reconciliation and P-Value Calculation
# Creator: Jorge Medina-Duran
# Date: Apr/27/2025
#
# Description:
# This script automates the execution of Empress reconciliations and p-value calculations
# across multiple subsampled parasite–host association matrices and cost schemes.
#
# Outputs:
# - CSV reconciliation files
# - SVG p-value plots
#
# Requirements:
# - Empress installed and accessible through specified PYTHON_EXEC path.
# - Prepared host/parasite trees and subsampled association matrices.
#############################################################

# Set up Python environment for Empress
PYTHON_EXEC="/home/user/miniconda3/envs/empress_env/bin/python"  # Change path for python accordingly

# Define paths
data_path="/mnt/c/YOUR_PATH/output" # Change path for python accordingly
subsampled_interactions="$data_path/subsampled_matrix"
results_path="$data_path/results_empress"

# Create results directory if it doesn't exist
mkdir -p "$results_path"

# Define paths to host and parasite trees
host_tree="$data_path/trimmed_tree/romalea_tree_trim.tre"
parasite_tree="$data_path/trimmed_tree/gregarine_tree_trim.tre"

# Define all cost schemes to test (d t l)
cost_schemes=(
  "1 1 1"
  "1 8 1"
  "2 1 2"
  "2 4 1"
  "4 1 1"
)

# Main loop over cost schemes
for scheme in "${cost_schemes[@]}"; do
  read d_rate t_rate l_rate <<< "$scheme"

  echo "------------------------------------------------------------"
  echo "Running cost scheme: d=$d_rate, t=$t_rate, l=$l_rate"
  echo "------------------------------------------------------------"

  # Loop over 50 subsampled association matrices
  for i in $(seq 1 50); do
    association_matrix="$subsampled_interactions/associations_${i}.mapping"
    reconciliation_file="$results_path/reconciliation_d${d_rate}_t${t_rate}_l${l_rate}_associations_${i}.csv"

    if [[ -f "$reconciliation_file" ]]; then
      echo "Skipping reconciliation for simulation $i (already exists)."
    else
      echo "Running reconciliation for simulation $i..."
      "$PYTHON_EXEC" empress_cli.py reconcile "$host_tree" "$parasite_tree" "$association_matrix" \
        -d $d_rate -t $t_rate -l $l_rate --csv "$reconciliation_file"
    fi
  done

  # Loop again for p-value calculations
  for i in $(seq 1 50); do
    pvalue_file="$results_path/res_d${d_rate}_t${t_rate}_l${l_rate}_associations_${i}.svg"

    if [[ -f "$pvalue_file" ]]; then
      echo "Skipping p-value calculation for simulation $i (already exists)."
    else
      echo "Running p-value calculation for simulation $i..."
      "$PYTHON_EXEC" empress_cli.py p-value "$host_tree" "$parasite_tree" "$subsampled_interactions/associations_${i}.mapping" \
        -d $d_rate -t $t_rate -l $l_rate --n 1000 --outfile "$pvalue_file"
    fi
  done

done
