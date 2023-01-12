#!/bin/bash
#SBATCH --job-name=gwas
#SBATCH --partition=jic-short,jic-medium,nbi-long,jic-long,RG-Cristobal-Uauy,nbi-medium
#SBATCH --nodes=1
#SBATCH --cpus=1
#SBATCH --mem 30G
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/quirozj/09_watseq/03_haplotypes/hapGWAS/scripts/log/gwas.%N.%j.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/quirozj/09_watseq/03_haplotypes/hapGWAS/scripts/log/gwas.%N.%j.err # STDERR
#SBATCH --array=0-1

i=$SLURM_ARRAY_TASK_ID

trait='SRA'
group='WatSeq_Pangenome_ABD'
dpm='slid_dmp'
genome='whole'
MAF=3
genotypes='watkins'
# chrom='chr2A'

# reference='sy_mattis'
declare -a references=(\
# "arinalrfor" \
# "chinese" \
# "jagger" \
"julius" \
# "lancer" \
# "landmark" \
# "mace" \
# "norin61" \
# "spelta" \
# "stanley" \
# "sy_mattis" \
)

reference=${references[$i]}

base_dir='/jic/scratch/groups/Cristobal-Uauy/quirozj/09_watseq/03_haplotypes/hapGWAS'
mtx_dir=$base_dir/${reference}/${group}/${dpm}
scripts_dir=$base_dir/scripts
out_dir=${base_dir}/01_results/${trait}
mkdir -p $out_dir

singularity exec ~/tmp/quirozc/python3.img python3 $scripts_dir/source/RunAssociation_GLM.py \
-i $mtx_dir/${genome}_genome_${reference}_${group}_${dpm}_hapgwas_${genotypes}.tsv.gz \
-hd $scripts_dir/metadata/watkins_accessions.txt \
-p $scripts_dir/phenotypes/${trait}.txt \
-pv 1 \
-c 0.2 \
-mc ${MAF} \
-dim 10 \
-pca $scripts_dir/metadata/watkins_pca.tsv \
-o $out_dir/${genome}_genome_${reference}_${trait}_${dpm}_gwas_${MAF}.tsv