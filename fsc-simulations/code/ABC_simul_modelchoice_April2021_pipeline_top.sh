#!/bin/bash -l
# NOTE the -l flag!
#SBATCH -J SCENARIO_START-TO
#SBATCH -o SCENARIO_START-TO.output
#SBATCH -e SCENARIO_START-TO.output
# Default in slurm
#SBATCH --mail-user user@email
#SBATCH --mail-type=FAIL,END
# Request 72 hours run time
#SBATCH -t 96:0:0
#SBATCH -A snic-project
#
#SBATCH -p core -n 2
# NOTE: You must not use more than 6GB of memory per core
# 1 core suffice for the migration with pulse simulations
# Time should be adapted to the number and nature of simulations. Pulse goes faster than continuous.

scenario=SCENARIO
start=START
end=TO

