#!/bin/csh

#Submit this script with: sbatch thefilename

#SBATCH --time=0:10:00   # walltime
##SBATCH --nodes=4    # number of nodes (40? cores per node)
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --mem=80G   # memory per CPU core
#SBATCH --error=master_it.e   # stderr file
#SBATCH --output=master_it.o   # stdout file
#SBATCH --exclusive
##SBATCH --hint=nomultithread

###***#SBATCH --job-name=nr02   # job name
###***#SBATCH --partition=   # queue to run in

echo "Starting" $SLURM_JOB_ID `date`
echo "Initiated on `hostname`"
echo ""
cd "$SLURM_SUBMIT_DIR"           # connect to working directory of sbatch

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
srun python calc_pyadjoint_misfit_obs.py
echo "Done" `date`
exit
