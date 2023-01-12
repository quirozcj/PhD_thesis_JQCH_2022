import io
import sys
from Phenotype_GLM import Phenotype
import pandas as pd
from sklearn.decomposition import PCA
import statsmodels.api as sm
import numpy as np
import math
from datetime import datetime
from bitarray import bitarray 

startTime = datetime.now()

class KmerProjection(object):
    def __init__(self, phenotype, haploMatrix, headerFile, snpFile, pcaDimensions, cor_threshold, pval_threshold, min_count, pca):
        self.associationMatrix = {}
        self.phenotype = phenotype
        self.cor_threshold = cor_threshold
        self.pval_threshold = pval_threshold
        self.min_count = min_count
        # self.presenceMatrix = presenceMatrix
        self.presenceMatrix = self.build_matrix(haploMatrix)
        # print(self.presenceMatrix)
        self.headerFile = headerFile
        self.snpFile = snpFile
        if pca is not None:
            self.pcaDF = pca
        else:
            self.pcaDF = self.computePCA(pcaDimensions)
        self.nullDF = self.pcaDF.copy()
        self.nullDF['bias'] = np.ones(self.nullDF.shape[0]).reshape(-1,1)
        self.null_results = sm.OLS(self.phenotype.phenoScores_series, self.nullDF).fit()

    def build_matrix(self, file):
        df = pd.read_csv(file, delimiter='\t')
        cols=list(df.columns[5:])
        df['presence'] = df[cols].apply(lambda row: ''.join(row.values.astype(str)), axis=1)
        matrix = df[['chr','start','end','haplotype','presence']]
        return matrix

    def computePCA(self, n_dimensions):
        snpDF = pd.read_csv(self.snpFile, sep = "\t", index_col=0)
        # retain only those accessions in SNP markers matrix for which phenotype scores are stored
        snpDF_reduced = snpDF.loc[self.phenotype.phenoScores_series.index]
        print("Computing PCA with " + str(n_dimensions) + " components.")
        pca = PCA(n_components = n_dimensions).fit_transform(snpDF_reduced)
        pcaDF= pd.DataFrame(pca, index = snpDF_reduced.index)
        pcaDF.columns = ['PCA_dim_' + str(col+1) for col in pcaDF.columns]
        # pcaDF.to_csv('tauchii_pca.tsv', sep = "\t" )
        return pcaDF

    def getCorrelation(self, presenceString, accessions_dict, pheno_sum, rho_pheno, n_accessions):
        presence_bitarray = bitarray(presenceString)
        presence_sum = 0.0
        presence_sqsum = 0.0
        presence_pheno_dot_product = 0.0
        for accession in self.phenotype.phenoScores_dict:
            index = accessions_dict[accession]
            if presence_bitarray[index]:
                presence_sum += 1
                presence_sqsum += 1
                presence_pheno_dot_product += self.phenotype.phenoScores_dict[accession]
        # if k-mer is present/absent in more than min_count, k-mer is given an correlation score of 0, ensuring that it is pre-filtered.
        if presence_sum < self.min_count or presence_sum > n_accessions-self.min_count:
            return 0
        rho_presence = math.sqrt(n_accessions*presence_sqsum - presence_sum*presence_sum)
        correlationScore = (n_accessions*presence_pheno_dot_product - presence_sum*pheno_sum)/(rho_presence*rho_pheno)    
        return correlationScore
   
    def getGLMpval(self, presenceString, header):
        presence_array = np.array(list(presenceString)).astype(int)
        presence_series = pd.Series(presence_array, index = header)
        pca_presenceDF = self.nullDF.copy()
        pca_presenceDF['presence'] = presence_series.loc[self.phenotype.phenoScores_series.index]
        # predict phenotype scores based on presence/absence of k-mer, with PCA dimensions as covariates        
        results = sm.OLS(self.phenotype.phenoScores_series, pca_presenceDF).fit()    
        # pvalue for nested models
        pval = results.compare_lr_test(self.null_results)[1]
        return -np.log10(pval)

    def readMatrix_GLM(self):
        accessions_dict = {}
        with open(self.headerFile, 'r') as h:
            header = [l.strip() for l in h]
        for i in range(len(header)):
            accessions_dict[header[i]] = i
        phenoScores_dict = self.phenotype.phenoScores_dict
        print("Correlation pre-filtering with a threshold of " + str(self.cor_threshold))
        # calculate phenotype parameters used in correlation calculation. reduces computation time (not by much).
        pheno_sum = 0.0
        pheno_sqsum = 0.0
        for accession in phenoScores_dict:
            pheno_sum += phenoScores_dict[accession]
            pheno_sqsum += phenoScores_dict[accession]*phenoScores_dict[accession]
        n_accessions = len(phenoScores_dict)
        rho_pheno = math.sqrt(n_accessions*pheno_sqsum - pheno_sum*pheno_sum)
        # rho_pheno is 0 if and only if phenotype scores are all equal. In that case, no association mapping is possible, therefore exit the program. 
        if rho_pheno == 0:
            print("Phenotype scores used in the correlation calculation cannot all be same.")
            sys.exit()    

        # file_in=pd.read_csv(self.presenceMatrix, delimiter='\t')
        file_in = self.presenceMatrix.copy()
        for index, row in file_in.iterrows():
            haplotype = row['haplotype']
            
            bool_haplotype = haplotype in self.associationMatrix
            if bool_haplotype:
                correlation = self.getCorrelation(row['presence'], accessions_dict, pheno_sum, rho_pheno, n_accessions)
                if abs(correlation) > self.cor_threshold:
                    pvalF = self.getGLMpval(row['presence'], header)
                    if pvalF > self.pval_threshold:
                        self.associationMatrix[haplotype] = [int(correlation*100),int(pvalF*100)]                             

    def readAssembly_GLM(self):
        # file_in=pd.read_csv(self.presenceMatrix, delimiter='\t')
        file_in = self.presenceMatrix.copy()
        for index, row in file_in.iterrows():
            haplotype = row['haplotype']
            self.associationMatrix[haplotype] = [0,-1]

    def writeAssociationScore_GLM(self, outputFile):
        out = open(outputFile, 'w')
        out.write('chromosome'+"\t"+'start'+"\t"+'end'+"\t"+'haplotype'+"\t"+'association'+"\t"+'correlation'+"\t"+'count'+"\n")
        
        # file_in=pd.read_csv(self.presenceMatrix, delimiter='\t')
        file_in = self.presenceMatrix.copy()
        for index, row in file_in.iterrows():
            chromosome = row['chr']
            start = row['start']
            end = row['end']
            haplotype = row['haplotype']
            
            h= {}
            correlation = self.associationMatrix[haplotype][0]/100
            associationScore  = self.associationMatrix[haplotype][1]/100
            if abs(correlation) > self.cor_threshold and associationScore > self.pval_threshold:
                if correlation < 0:
                    associationScore = -1*associationScore
                num = 0
                if associationScore in h:
                    num = h[associationScore][0]
                num += 1
                h[associationScore] = [correlation, num]
            for associationScore, ele in h.items():
                out.write(chromosome+"\t"+str(start)+"\t"+str(end)+"\t"+haplotype+"\t"+str(abs(associationScore))+"\t"+str(ele[0])+"\t"+str(ele[1])+"\n")
        out.close()
        print(datetime.now() - startTime)