#! /bin/bash -l

bamfile_dir="/athena/angsd/scratch/yus4008/project/dataset/trim_alignments/"
gtffile="/athena/angsd/scratch/yus4008/project/dataset/hg38_genome/hg38.ncbiRefSeq.gtf"
resdir="/athena/angsd/scratch/yus4008/project/dataset/qc_analysis/qorts/results/"

mamba activate qorts

for response in uninfected/ symptomatic/
do
    bam_use_dir=${bamfile_dir}${response}
    for file in ${bam_use_dir}*.bam
    do
        acc=`echo $file | egrep -o "SRR([0-9]+)"`
        outdir="${resdir}${acc}/"
        if [ ! -d "$outdir" ]
        then
            mkdir $outdir
        fi
        qorts -Xmx16G QC --generatePlots --maxPhredScore 45 $file $gtffile $outdir
    done
done

mamba deactivate