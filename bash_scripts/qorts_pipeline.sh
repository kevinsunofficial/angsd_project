#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=qorts_pipeliine
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --mail-user=yus4008@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --output=qorts.out

bash -l make_decoder.sh
bash -l run_qorts.sh

mamba activate qorts
Rscript plot.R
mamba deactivate
