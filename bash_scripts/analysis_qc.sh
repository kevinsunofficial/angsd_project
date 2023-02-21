#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=align_seqs
#SBATCH --time=08:00:00
#SBATCH --mem=32G
#SBATCH --mail-user=yus4008@med.cornell.edu
#SBATCH --mail-type=ALL


script_dir="/home/yus4008/cmpb5004/project/angsd_project/bash_scripts/"
dataset_dir="/athena/angsd/scratch/yus4008/project/dataset/"
aln_dir="/athena/angsd/scratch/yus4008/project/dataset/alignments/"
qc_dir="${dataset_dir}analysis_qc/"
fastqc_dir="${qc_dir}fastqc/"
bamqc_dir="${qc_dir}bamqc/"
multiqc_dir="${qc_dir}multiqc/"

aligned=true

mamba activate angsd

for response in uninfected/ symptomatic/
do
    fastqc_use_dir="${dataset_dir}${response}"
    fastqc_out_dir="${fastqc_dir}${response}"
    for file in $fastqc_use_dir*
    do
        acc=`echo $file | egrep -o "SRR([0-9]+)_(1|2)"`
        echo -e "fastqc ${file} --noextract --outdir ${fastqc_out_dir}"
        fastqc $file --noextract --outdir $fastqc_out_dir
    done

    if [ "${aligned}" = true ]
    then
        bamqc_use_dir="${aln_dir}${response}"
        bamqc_out_dir="${bamqc_dir}${response}"
        echo -e "/softlib/apps/EL7/BamQC/bin/bamqc -o ${bamqc_out_dir} --noextract ${bamqc_use_dir}*.bam"
        /softlib/apps/EL7/BamQC/bin/bamqc -o $bamqc_out_dir --noextract $bamqc_use_dir*.bam
    fi
done

mamba deactivate


mamba activate multiqc

echo -e "multiqc -o ${multiqc_dir} ${qc_dir}"
multiqc -o $multiqc_dir $qc_dir

mamba deactivate
