#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=bamqc
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --mail-user=yus4008@med.cornell.edu
#SBATCH --mail-type=ALL


script_dir="/home/yus4008/cmpb5004/project/angsd_project/bash_scripts/"
dataset_dir="/athena/angsd/scratch/yus4008/project/dataset/"
aln_dir="/athena/angsd/scratch/yus4008/project/dataset/alignments/"
qc_dir="${dataset_dir}qc_analysis/"
bamqc_dir="${qc_dir}bamqc/"
multiqc_dir="${qc_dir}multiqc/bamqc/"


mamba activate angsd

for response in uninfected/ symptomatic/
do
    bamqc_use_dir="${aln_dir}${response}"
    bamqc_out_dir="${bamqc_dir}${response}"
    echo -e "/softlib/apps/EL7/BamQC/bin/bamqc -o ${bamqc_out_dir} --noextract ${bamqc_use_dir}*.bam"
    /softlib/apps/EL7/BamQC/bin/bamqc -o $bamqc_out_dir --noextract $bamqc_use_dir*.bam
done

mamba deactivate

mamba activate multiqc

echo -e "multiqc -o ${multiqc_dir} ${bamqc_dir}"
multiqc -o $multiqc_dir $bamqc_dir

mamba deactivate
