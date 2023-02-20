#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=align_seqs
#SBATCH --time=08:00:00
#SBATCH --mem=100G
#SBATCH --mail-user=yus4008@med.cornell.edu
#SBATCH --mail-type=ALL


script_dir="/home/yus4008/cmpb5004/project/angsd_project/bash_scripts/"
aln_dir="/athena/angsd/scratch/yus4008/project/dataset/alignments/"
qc_dir="${aln_dir}analysis_qc/"
bamqc_dir="${qc_dir}bamqc/"
multiqc_dir="${qc_dir}multiqc/"


mamba activate angsd

/softlib/apps/EL7/BamQC/bin/bamqc -o $bamqc_dir --noextract $aln_dir/*/*.bam

mamba deactivate

mamba activate multiqc

multiqc -o $multiqc_dir $aln_dir

mamba deactivate
