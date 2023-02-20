#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=20230212_align
#SBATCH --time=08:00:00
#SBATCH --mem=100G
#SBATCH --mail-user=yus4008@med.cornell.edu
#SBATCH --mail-type=ALL


script_dir="/home/yus4008/cmpb5004/project/angsd_project/bash_scripts/"
ref_dir="${dataset_dir}hg38_STARindex"
dataset_dir="/athena/angsd/scratch/yus4008/project/dataset/"
alignment_dir="${dataset_dir}alignments/"

limit=1

mamba activate angsd

for response in uninfected/ symptomatic/
do
    use_dir="${dataset_dir}${response}"
    current=0
    for file in "${use_dir}*"
    do
        if [[ $current < $limit ]]
        then
            acc=`echo $file | egrep -o "SRR([0-9]+)_(1|2)"`
            out_dir="${alignment_dir}${response}${acc}."
            out_file="${out_dir}Aligned.sortedByCoord.out.bam"
            echo $file $acc $out_dir $out_file
            if [ ! -f "${out_file}" ]
            then
                echo $out_file exists
                # STAR --runMode alignReads \
                #      --runThreadN 1 \
                #      --genomeDir $ref_dir \
                #      --readFilesIn $file \
                #      --readFilesCommand zcat \
                #      --outFileNamePrefix $out_dir \
                #      --outSAMtype BAM SortedByCoordinate
            fi
            # samtools index $out_file
        fi
    done
    current=$((current+1))
done

mamba deactivate
