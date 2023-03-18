#! /bin/bash -l

#SBATCH --partition=angsd_class
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=feature_counts
#SBATCH --time=48:00:00
#SBATCH --mem=64G
#SBATCH --mail-user=yus4008@med.cornell.edu
#SBATCH --mail-type=ALL
#SBATCH --output=feature_counts.out

bam_dir="/athena/angsd/scratch/yus4008/project/dataset/trim_alignments/"
gtf="/athena/angsd/scratch/yus4008/project/dataset/hg38_genome/hg38.ncbiRefSeq.gtf"
fc_dir="/athena/angsd/scratch/yus4008/project/dataset/featureCounts/"

mamba activate angsd

featureCounts -a $gtf -t exon -g gene_id -o $fc_dir $bam_dir*/*.bam

mamba deactivate