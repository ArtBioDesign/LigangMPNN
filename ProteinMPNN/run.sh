#!/bin/bash

# Function to print usage


echo "Received arguments: $@"


site_pdb=$1
complex_pdb=$2
out_dir=$3

# Define the environment and script paths
env="/opt/conda/envs/ligandmpnn_env"
ligandpnn="/workspace/pmpnn/LigandMPNN/run.py"


# Print received arguments
echo "Received arguments:"
echo "site_pdb: $site_pdb"
echo "complex_pdb: $complex_pdb"
echo "out_dir: $out_dir"

# Run Python script
/opt/conda/envs/ligandmpnn_env/bin/python   /workspace/pmpnn/run.py --site "$site_pdb" --complex_pdb "$complex_pdb" --out_dir "$out_dir" --env "$env" --ligandpnn "$ligandpnn"

