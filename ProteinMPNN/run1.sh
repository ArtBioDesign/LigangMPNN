#!/bin/bash

# Function to print usage
usage() {
    echo "Usage: $0"
    exit 1
}

# Function to check if required files and directories exist
check_requirements() {
    local required_files=("$site_pdb" "$complex_pdb")
    local required_dirs=("$out_dir")

    for file in "${required_files[@]}"; do
        if [ ! -f "$file" ]; then
            echo "Error: Required file $file does not exist."
            exit 1
        fi
    done

    for dir in "${required_dirs[@]}"; do
        if [ ! -d "$dir" ]; then
            echo "Error: Required directory $dir does not exist."
            exit 1
        fi
    done
}

# Set environment variables
set_env_vars() {
    export TMPDIR=/tmp
    export PATH="/opt/conda/envs/ligandmpnn_env/bin:$PATH"
}

# Print received arguments
print_args() {
    echo "Received arguments: $@"
}

# Define paths
BASE_DIR=$(pwd)
site_pdb="/workspace/pmpnn/inputs/site.pdb"
complex_pdb="/workspace/pmpnn/inputs/complex.pdb"
out_dir="/tmp/results"
env="/opt/conda/envs/ligandmpnn_env"
ligandpnn="/workspace/pmpnn/LigandMPNN/run.py"

# Check requirements
check_requirements

# Set environment variables
set_env_vars

# Print arguments
print_args "$@"

# Run Python script
python /workspace/pmpnn/run.py --site_pdb "$site_pdb" --complex_pdb "$complex_pdb" --out_dir "$out_dir" --env "$env" --ligandpnn "$ligandpnn"

