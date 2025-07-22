#! /bin/bash

#SBATCH --partition=Orion
#SBATCH --job-name=my-first-cluster-run
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=0:02:00

echo "======================================================"
echo "Start Time  : $(date)"
echo "Submit Dir  : $SLURM_SUBMIT_DIR"
echo "Job ID/Name : $SLURM_JOBID / $SLURM_JOB_NAME"
echo "Node List   : $SLURM_JOB_NODELIST"
echo "Num Tasks   : $SLURM_NTASKS total [$SLURM_NNODES nodes @ $SLURM_CPUS_ON_NODE CPUs/node]"
echo "Num GPUs    : $NUM_GPUS $GPU_MDL"
echo "======================================================"
echo ""

python3 ICU-preparedness-cluster.py

echo ""
echo "======================================================"
echo "End Time   : $(date)"
echo "======================================================"

