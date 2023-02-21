#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=create_index
#SBATCH --time=08:00:00
#SBATCH --mem=100G
#SBATCH --mail-user=yus4008@med.cornell.edu
#SBATCH --mail-type=ALL


script_dir="/home/yus4008/cmpb5004/project/angsd_project/bash_scripts/"
dataset_dir="/athena/angsd/scratch/yus4008/project/dataset/"
genome_dir="${dataset_dir}hg38_genome/"
genome_seq="${genome_dir}hg38.fa"
genome_annot="${genome_dir}hg38.ncbiRefSeq.gtf"
genome_seq_gz="${genome_seq}.gz"
genome_annot_gz="${genome_annot}.gz"
ref_dir="${dataset_dir}hg38_STARindex"


if [ ! -f $genome_seq ]
then
    if [ ! -f $genome_seq_gz ]
    then
        echo -e "wget -O ${genome_seq_gz} http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
        wget -O $genome_seq_gz "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
    fi
    echo -e "gzip -d ${genome_seq_gz}"
    gzip -d $genome_seq_gz
fi

if [ ! -f $genome_annot ]
then
    if [ ! -f $genome_annot_gz ]
    then
        echo -e "wget -O ${genome_seq_gz} http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz"
        wget -O $genome_annot_gz "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz"
    fi
    echo -e "gzip -d ${genome_annot_gz}"
    gzip -d $genome_annot_gz
fi


mamba activate angsd

echo "mamba activated, switching directory to: ${ref_dir}"

cd $ref_dir

echo -e "STAR \t --runMode genomeGenerate \n
            \t\t --runThreadN 1 \n
            \t\t --genomeDir ${ref_dir} \n
            \t\t --genomeFastaFiles ${genome_seq} \n
            \t\t --sjdbGTFfile ${genome_annot} \n
            \t\t --sjdbOverhang 99"

STAR --runMode genomeGenerate \
     --runThreadN 1 \
     --genomeDir $ref_dir \
     --genomeFastaFiles $genome_seq \
     --sjdbGTFfile $genome_annot \
     --sjdbOverhang 99

cd $script_dir

echo "switching directory to: ${script_dir}, deactivate mamba"

mamba deactivate
