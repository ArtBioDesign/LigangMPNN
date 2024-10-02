#!/bin/bash

#SBATCH --job-name=pmpnn
#SBATCH --partition=qgpu_3090
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem=50G
#SBATCH --gres=gpu:1
#SBATCH --mail-type=end
#SBATCH --mail-user=YOU@MAIL.COM
#SBATCH --output=%j.out
#SBATCH --error=%j.err




singularity exec --nv  --writable-tmpfs pmpnn_latest.sif  /opt/conda/envs/ligandmpnn_env/bin/python /hpcfs/fhome/yangchh/pmpnn/ProteinMPNN/run.py --site /hpcfs/fhome/yangchh/pmpnn/inputs/1-8A1.csv --complex_pdb /hpcfs/fhome/yangchh/pmpnn/inputs/1-4CL-1.pdb --out_dir /hpcfs/fhome/yangchh/pmpnn/output --env /opt/conda/envs/ligandmpnn_env  --ligandpnn /hpcfs/fhome/yangchh/pmpnn/ProteinMPNN/LigandMPNN/run.py


