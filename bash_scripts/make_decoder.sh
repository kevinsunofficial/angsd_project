#! /bin/bash -l

bamfile_dir="/athena/angsd/scratch/yus4008/project/dataset/trim_alignments/"
dcdtxt="/athena/angsd/scratch/yus4008/project/dataset/qc_analysis/qorts/decoder.txt"

echo -e "unique.ID\tgroup.ID" > ${dcdtxt}

if [ -f "$dcdtxt" ]
then
    rm $dcdtxt
fi

for response in uninfected symptomatic
do
    bam_use_dir=${bamfile_dir}${response}/
    for file in ${bam_use_dir}*.bam
    do
        acc=`echo $file | egrep -o "SRR([0-9]+)"`
        echo -e "${acc}\t${response}" >> ${dcdtxt}
    done
done
