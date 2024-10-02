# LigangMPNN
## **Project Introduction**  
**LigangMPNN** is a protein design tool based on complexes, developed by the Baker group in 2023.

## Installation & Execution for HPC Deployment

### **1. Build the Environment**
```bash
conda env create -f environment.yml
```

### **2. Build image from docker file**
```shell
docker build -f dockerfile_env -t pmpnn_env:latest .
docker build -f dockerfile -t pmpnn:latest .
```

### **3. Save the Docker Image**
```shell
docker save -o pmpnn_latest.tar pmpnn:latest
```

### **4. Convert Docker Image to Singularity**
```shell
singularity build pmpnn_latest.sif docker-archive://pmpnn_latest.tar
```

### **5. Submit the Job to SLURM**
```shell
sbatch runhpc.sh
```