#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=20230212_align
#SBATCH --time=05:00:00
#SBATCH --mem=32G
#SBATCH --mail-user=yus4008@med.cornell.edu
#SBATCH --mail-type=ALL


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
        wget -O $genome_seq_gz "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
    fi
    gzip -d $genome_seq_gz
fi

if [ ! -f $genome_annot ]
then
    if [ ! -f $genome_annot_gz ]
    then
        wget -O $genome_annot_gz "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz"
    fi
    gzip -d $genome_annot_gz
fi


mamba activate angsd

echo "mamba activated"

STAR --runMode genomeGenerate --runThreadN 1 --genomeDir $ref_dir --genomeFastaFiles $genome_seq --sjdbGTFfile $genome_annot --sjdbOverhang 99

mamba deactivate