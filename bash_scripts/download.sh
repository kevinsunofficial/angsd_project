#! /bin/bash

entry_dir="../entry/"
ena_report="${entry_dir}filereport_read_run_PRJNA744408_tsv.txt"
sra_table="${entry_dir}SraRunTable.csv"
srr_list="${entry_dir}SRR_Acc_List.txt"

dataset_dir="/athena/angsd/scratch/yus4008/project/dataset/"
uninfected_dir="${dataset_dir}uninfected/"
symptomatic_dir="${dataset_dir}symptomatic/"

limit=6

while read acc
do
    response=`awk -v pat="${acc}" -F "," '{ if ($1 ~ pat && $2 <= 50 && $2 >= 30) print $4; }' $sra_table`
    if [ ! -z "$response" ]
    then
        response_dir="${dataset_dir}${response}/"
        current=`ls $response_dir | wc -l`
        if [[ $current -lt $limit ]]
        then
            fastq1="${response_dir}${acc}_1.fastq.gz"
            fastq2="${response_dir}${acc}_2.fastq.gz"
            line1=`egrep -o "ftp.sra.ebi.ac.uk/vol1/fastq/SRR150/[0-9]{3}/${acc}/${acc}_1.fastq.gz" $ena_report`
            line2=`egrep -o "ftp.sra.ebi.ac.uk/vol1/fastq/SRR150/[0-9]{3}/${acc}/${acc}_2.fastq.gz" $ena_report`
            
            if [ ! -z "$line1" ]
            then
                wget -O $fastq1 ftp://$line1
            fi
            if [ ! -z "$line2" ]
            then
                wget -O $fastq2 ftp://$line2
            fi
        fi
    fi
done < $srr_list
