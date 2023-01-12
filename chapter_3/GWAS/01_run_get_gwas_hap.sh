#!/bin/bash
#SBATCH --partition=jic-short,jic-medium,nbi-long,jic-long,RG-Cristobal-Uauy,nbi-medium
#SBATCH --nodes=1
#SBATCH --cpus=1
#SBATCH --mem 10G
#SBATCH -o /jic/scratch/groups/Cristobal-Uauy/quirozj/09_watseq/03_haplotypes/hapGWAS/scripts/log/gws_f.%N.%j.out # STDOUT
#SBATCH -e /jic/scratch/groups/Cristobal-Uauy/quirozj/09_watseq/03_haplotypes/hapGWAS/scripts/log/gws_f.%N.%j.err # STDERR
#SBATCH --job-name=gws_f
#SBATCH --array=0-20

i=$SLURM_ARRAY_TASK_ID

declare -a chromosomes=()
for n in {1..7..1}; do
	for l in A B D; do
		chromosomes+=(chr${n}${l});
	done;
done

# chromosome='chr7B'
# declare -a chromosomes=('chr1A','chr7B')
chr_i=$(($i%21))
chromosome=${chromosomes[$chr_i]}

# chromosome=${chromosomes[$i]}
# reference=sy_mattis
declare -a references=(\
"arinalrfor" \
# "chinese" \
# "jagger" \
# "julius" \
# "lancer" \
# "landmark" \
# "mace" \
# "norin61" \
# "spelta" \
# "stanley" \
# "sy_mattis" \
)

ref_i=$(($i/21))
reference=${references[$ref_i]}

group='WatSeq_Pangenome_ABD'
dmps='gwas_dmp'
# format='plink'
format='hapgwas'
genome='whole'
genotypes='watkins_addlines'

out_dir=/jic/scratch/groups/Cristobal-Uauy/quirozj/09_watseq/03_haplotypes/hapGWAS/${reference}/${group}/${dmps}
mkdir -p ${out_dir}

singularity exec ~/tmp/quirozc/python3.img python3 get_gwas_hap.py \
-f ${group}_${dmps}.tsv \
-r ${reference} \
-c ${chromosome} \
-l ${genotypes}.tsv \
-p chr_prefix.tsv \
-o ${out_dir}/${chromosome}_${reference}_${group}_${dmps}_${format}.tsv.gz
# -s ${format} \ # note, for format of GAPIT

wait
# combine individual
cd ${out_dir}
zcat *${format}.tsv.gz |awk '!seen[$0]++'|tr ' ' '\t' | gzip > ${genome}_genome_${reference}_${group}_${dmps}_${format}_${genotypes}.tsv.gz
# zcat *watkins.tsv.gz |awk '!seen[$0]++'|tr ' ' '\t' | gzip > whole_genome_${reference}_${group}_${dmps}_hapgwas_watkins_${format}_combined.tsv.gz

# zcat *${genome}_${reference}* |awk '!seen[$0]++'|tr ' ' '\t' | gzip > ${genome}_genome_${reference}_${group}_${dmps}_watkins_${format}.tsv.gz