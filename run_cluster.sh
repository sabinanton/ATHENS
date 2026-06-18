#!/bin/sh
# 
#SBATCH --job-name="sw_he_02"
#SBATCH --partition=compute-p2
#SBATCH --time 120:00:00
#SBATCH --ntasks=250
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --account=research-ae-spe

module load 2024r1
module load openmpi

SIM_ID=""

# Parse command-line arguments
for i in "$@"
do
case $i in
    -id=*)
    SIM_ID="${i#*=}"
    shift # past argument=value
    ;;
esac
done

srun ./athens --batch --project=Swarm_LookUp_Helium_02 --id=$SIM_ID