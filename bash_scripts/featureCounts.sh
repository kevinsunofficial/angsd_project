#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=feature_counts
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --mail-user=yus4008@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --output=feature_counts.out

bam_dir="/athena/angsd/scratch/yus4008/project/dataset/trim_alignments/"
gtf="/athena/angsd/scratch/yus4008/project/dataset/hg38_genome/hg38.ncbiRefSeq.gtf"
fc_dir="/athena/angsd/scratch/yus4008/project/dataset/featureCounts/feature_counts.txt"

mamba activate angsd

featureCounts -p --countReadPairs -a $gtf -t exon -g gene_id -o $fc_dir $bam_dir*/*.bam

mamba deactivate
