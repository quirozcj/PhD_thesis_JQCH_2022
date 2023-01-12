from KmerProjection_GLM import KmerProjection
from Phenotype_GLM import Phenotype
import sys 
import argparse
import os.path
from pathlib import Path
import pandas as pd

if __name__ == '__main__':
    if sys.version_info[0] < 3:
        sys.stderr.write('Python version 3 or above required') 
    parser = argparse.ArgumentParser(description = "Kmer Association mapping with Generalized Linear Model")
    inputMatrix = parser.add_argument_group('Kmer presence/absence matrix')
    inputMatrix.add_argument('-i', '--inputmatrix',  required=True, help='presence/absence matrix of haplotypes')
    inputMatrix.add_argument('-hd', '--header',  required=True, help='file header for the presence/absence matrix of kmers')
    phenotype = parser.add_argument_group('Phenotype')
    phenotype.add_argument('-p','--phenotype', required=True, help='Path to phenotype file')
    SNPmarker = parser.add_argument_group('SNP markers matrix')
    SNPmarker.add_argument('-s','--snp', default = None, help='SNP markers to compute PCA and correct for population structure')
    SNPmarker.add_argument("-dim", "--pcadimensions", type=int, default = 3, help="Number of significant PCA dimensions retained for regression analysis")
    output = parser.add_argument_group('Output')
    output.add_argument('-o','--output', required=True, help='ouput file')
    parser.add_argument("-u", "--usable", default = None, help="usable accessions")
    parser.add_argument("-c", "--correlationthreshold", type=float, default = 0.2, help="haplotypes whose correlation with phenotype is greater than this value are retained for regression")
    parser.add_argument("-mc", "--mincount", type=int, default = 4, help="haplotypes retained for regression analysis which are present/absent in more than this number of accessions.")
    parser.add_argument("-pv", "--pvalthreshold", type=float, default = 3.0, help="only those k-mers with log10 of p-value greater than this value are retained")
    parser.add_argument("-st", "--stackman", default = None, help="convert phenotype scores from Stackman's IT to AgRenSeq scores", action="store_true")
    parser.add_argument("-per", "--permute", help="permute phenotype", action="store_true")
    parser.add_argument("-sub", "--subsample", type=int, help="random subsample of this size")
    parser.add_argument("-pca", "--pca_file", default = None)
    
    
    args = parser.parse_args()
    inputMatrix_f = args.inputmatrix
    header_f = args.header
    phenotype_filename = args.phenotype
    stackman = args.stackman

    try:
        phenotype_f = Phenotype(phenotype_filename, args.stackman)
    except FileNotFoundError:
        print("\nFile " + phenotype_filename + " does not exist.")
        sys.exit()
    if args.usable is not None:
        try:
            phenotype_f.selectAccessions(args.usable)
        except FileNotFoundError:
            print("\nFile " + args.usable + " does not exist.")
            sys.exit()
    if args.subsample is not None:       
        phenotype_f.selectRandomAccessions(args.subsample)
    if args.permute:
        phenotype_f.permutePhenotype()  
    print("\nThe number of accessions used for association are: " + str(len(phenotype_f.phenoScores_dict)))

    snp_f = args.snp        
    pca_dimensions = args.pcadimensions
    cor_threshold = args.correlationthreshold
    min_count = args.mincount
    pval_threshold = args.pvalthreshold
    output_filename = args.output
    # if os.path.isfile(output_filename):
    #     print('\nFile ' + output_filename + ' already exists')

    accessions = list(pd.read_csv(phenotype_filename, sep = "\t", header=None)[0])
    pca_file = args.pca_file
    # filter_pca = args.pca_file
    if pca_file is not None:    
        pca_file = pd.read_csv(pca_file, sep = "\t", index_col=0)
        filter_pca = pca_file.T.filter(items=accessions, axis=1).T
    projection = KmerProjection(phenotype_f, inputMatrix_f, header_f, snp_f, pca_dimensions, cor_threshold, pval_threshold, min_count, filter_pca)
			
    projection.readAssembly_GLM()
    projection.readMatrix_GLM()
    projection.writeAssociationScore_GLM( output_filename )