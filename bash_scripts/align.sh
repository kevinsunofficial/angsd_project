#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=align_seqs
#SBATCH --time=08:00:00
#SBATCH --mem=64G
#SBATCH --mail-user=yus4008@med.cornell.edu
#SBATCH --mail-type=ALL


script_dir="/home/yus4008/cmpb5004/project/angsd_project/bash_scripts/"
dataset_dir="/athena/angsd/scratch/yus4008/project/dataset/"
ref_dir="${dataset_dir}hg38_STARindex"
alignment_dir="${dataset_dir}trim_alignments/"

limit=1

mamba activate angsd

for response in uninfected/ symptomatic/
do
    use_dir="${dataset_dir}${response}trim/"
    current=0
    for file in $use_dir*_1_val_1.fq.gz
    do
        if [[ $current -lt $limit ]]
        then
            acc=`echo $file | egrep -o "SRR([0-9]+)"`
            file2="${use_dir}${acc}_2_val_2.fq.gz"
            out_dir="${alignment_dir}${response}${acc}."
            out_file="${out_dir}Aligned.sortedByCoord.out.bam"
            out_flagstat="${out_dir}flagstats"
            out_stats="${out_dir}stats"
            if [ ! -f "${out_file}" ]
            then
                echo -e "STAR --runMode alignReads --readFilesIn $acc (${current})"
                STAR --runMode alignReads \
                     --runThreadN 1 \
                     --genomeDir $ref_dir \
                     --readFilesIn $file $file2 \
                     --readFilesCommand zcat \
                     --outFileNamePrefix $out_dir \
                     --outSAMtype BAM SortedByCoordinate
                echo -e "samtools index ${out_file}"
                samtools index $out_file
            fi
            if [ ! -f "${out_flagstat}" ]
            then
                echo -e "samtools flagstat ${out_file} > ${out_flagstat}"
                samtools flagstat $out_file > $out_flagstat
            fi
            if [ ! -f "${out_stats}" ]
            then
                echo -e "samtools stats ${out_file} > ${out_stats}"
                samtools stats $out_file > $out_stats
            fi
            current=$((current+1))
        fi
    done
done

mamba deactivate
