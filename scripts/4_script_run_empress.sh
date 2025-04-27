#!/bin/bash

# Set up Python environment for Empress
#export PYTHONUSERBASE=$WORK/cophylo_signal/emPress/script/.local_py3.7.10
#export PATH=$PYTHONUSERBASE/bin:$PATH
PYTHON_EXEC=/usr/bin/python3.7


# Define paths
data_path="/mnt/c/Users/beto2/Documents/cophylogeny_project/output" # change path to your data directory
subsampled_interactions="$data_path/subsampled_matrix"
results_path="$data_path/results_empress"
mkdir -p "$results_path"

# Set reconciliation parameters
d_rate=4
t_rate=1
l_rate=1

# Define host and parasite phylogenies (unchanged across runs)
host_tree="$data_path/trimmed_tree/romalea_tree_trim.tre"
parasite_tree="$data_path/trimmed_tree/gregarine_tree_trim.tre"

# Iterate over 50 association matrices
for i in $(seq 1 50)
do
    association_matrix="$data_path/associations_${i}.mapping"

    # Reconciliation file name
    reconciliation_file="$results_path/reconciliation_d${d_rate}_t${t_rate}_l${l_rate}_associations_${i}.csv"

    # Run Empress Reconciliation if output does not exist
    if test -f "$reconciliation_file"; then
        echo "Skipping reconciliation for simulation ${i}, already exists."
    else
        echo "Running reconciliation for simulation ${i}"
        $PYTHON_EXEC empress_cli.py reconcile "$host_tree" "$parasite_tree" "$association_matrix" \
        -d $d_rate -t $t_rate -l $l_rate --csv "$reconciliation_file"
    fi
done

# Run Empress P-Value Calculation for Each Association Matrix
for i in $(seq 1 50)
do
    pvalue_file="$results_path/res_d${d_rate}_t${t_rate}_l${l_rate}_associations_${i}.svg"

    if test -f "$pvalue_file"; then
        echo "Skipping p-value calculation for simulation ${i}, already exists."
    else
        echo "Running p-value calculation for simulation ${i}"
        $PYTHON_EXEC empress_cli.py p-value "$host_tree" "$parasite_tree" "$data_path/associations_${i}.mapping" \
        -d $d_rate -t $t_rate -l $l_rate --n 1000 --outfile "$pvalue_file"
    fi
done
