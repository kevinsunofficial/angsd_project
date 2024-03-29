#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=fastqc_trim
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --mail-user=yus4008@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --output=fastqc_trim.out


script_dir="/home/yus4008/cmpb5004/project/angsd_project/bash_scripts/"
dataset_dir="/athena/angsd/scratch/yus4008/project/dataset/"
aln_dir="/athena/angsd/scratch/yus4008/project/dataset/alignments/"
qc_dir="${dataset_dir}qc_analysis/"
fastqc_dir="${qc_dir}fastqc/"
untrim_fastqc_dir="${fastqc_dir}untrim/"
trim_fastqc_dir="${fastqc_dir}trim/"
multiqc_dir="${qc_dir}multiqc/"
untrim_multiqc_dir="${multiqc_dir}fastqc_untrimmed"
trim_multiqc_dir="${multiqc_dir}fastqc_trimmed"

mamba activate angsd

for response in uninfected/ symptomatic/
do
    # echo -e "Perform fastqc on raw reads"
    fastqc_use_dir="${dataset_dir}${response}"
    fastqc_out_dir="${untrim_fastqc_dir}${response}"
    # for file in $fastqc_use_dir*.fastq.gz
    # do
    #     acc=`echo $file | egrep -o "SRR([0-9]+)_(1|2)"`
    #     echo -e "fastqc ${file} --noextract --outdir ${fastqc_out_dir}"
    #     fastqc $file --noextract --outdir $fastqc_out_dir
    # done

    echo -e "Perform trim-galore"
    mamba activate trim-galore
    trim_out_dir="${fastqc_use_dir}trim/"
    for file in $fastqc_use_dir*_1.fastq.gz
    do
        acc=`echo $file | egrep -o "SRR([0-9]+)"`
        file2="${fastqc_use_dir}${acc}_2.fastq.gz"
        echo -e "trim_galore --illumina ${file} ${file2} --output_dir ${trim_out_dir} --stringency 13"
        trim_galore --illumina --paired --output_dir $trim_out_dir --stringency 13 $file $file2
    done
    mamba deactivate

    echo -e "Perform fastqc on trimmed reads"
    fastqc_use_dir="${dataset_dir}${response}trim/"
    fastqc_out_dir="${trim_fastqc_dir}${response}"
    for file in $fastqc_use_dir*.fq.gz
    do
        acc=`echo $file | egrep -o "SRR([0-9]+)_(1|2)"`
        echo -e "fastqc ${file} --noextract --outdir ${fastqc_out_dir}"
        fastqc $file --noextract --outdir $fastqc_out_dir
    done
done

mamba deactivate

mamba activate multiqc

# echo -e "multiqc -o ${untrim_multiqc_dir} ${untrim_fastqc_dir}"
# multiqc -o $untrim_multiqc_dir $untrim_fastqc_dir
echo -e "multiqc -o ${trim_multiqc_dir} ${trim_fastqc_dir}"
multiqc -o $trim_multiqc_dir $trim_fastqc_dir

mamba deactivate
