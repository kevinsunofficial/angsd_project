#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=multiqc
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --mail-user=yus4008@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --output=multiqc.out


multiqc_dir="/athena/angsd/scratch/yus4008/project/dataset/"
alndir="${multiqc_dir}trim_alignments/"
fcdir="${multiqc_dir}featureCounts/"
qcdir="${multiqc_dir}qc_analysis/"
resdir="/athena/angsd/scratch/yus4008/project/dataset/qc_analysis/multiqc/"

mamba activate multiqc

echo -e "multiqc ${alndir} ${fcdir} ${qcdir} -o ${resdir}"
multiqc $alndir $fcdir $qcdir -o $resdir

mamba deactivate