# Aim: define the k-mer size for haplotype bilding
#### NOTE: in another experiment I'll test a different k-mer counter:
#### https://academic.oup.com/bioinformatics/article/33/17/2759/3796399
#### KMC 3: counting and manipulating k-mer statistics

## K-mer size to test:
```
k=11
k=21
k=31
k=41
k=51
k=61
k=71
k=81
k=91
k=101
```
Working directory:
```
/jic/scratch/groups/Cristobal-Uauy/quirozj/00_jellies/06_define_kmer_size_genome
# in my local computer:
/Users/quirozc/Documents/01_hap_temporary/define_kmer_size

# script:
jellyfish_define_kmer_size.sh
# comand in the script:

NOTE: do not use -C when genome assemblies, since canonical kmers (reverse complement will be counted as a the same)

jellyfish count -C -m 11 -s 3G --disk --out-counter-len 1 --disk -o $outdir/kmer_size/chinese_spring_k11.jf -t 5 $indir/161010_Chinese_Spring_v1.0_pseudomolecules

```

The time requered for running the previos command for each k-mer and the chinese spring genom size is usin 5 thereads/CPUs:

```
# k=11
# k=21 # for this size takes ~ 30 min for chinese spring assembly using 10 threads
# k=31 # 1 hour
# k=41 # 1 1/2 hour
# k=51
# k=61 # 2 hours
# k=71 # 2 1/2 hours
# k=81 # 3 hours
# k=91 # 3 1/2 
# k=101
# k=111
# k=121
# k=131
# k=141
# k=151
```

### get the histogram using jellysfish:
```
jellyfish histo chinese_spring_k31.jf -t 32 -f -o chinese_spring_k31.histo.txt

# get the stats:

jellyfish stats chinese_spring_k31.jf > chinese_spring_k31.stats.txt

```


scripts and further analysis are in:

```sh
/Volumes/quirozj/09_watseq/14_kmer_histo
```
output of the files histo combined are:
```sh
/Volumes/quirozj/09_watseq/14_kmer_histo/histo
```











