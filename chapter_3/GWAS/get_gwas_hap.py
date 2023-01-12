import argparse
import pandas as pd
import numpy as np
from warnings import filterwarnings
filterwarnings('ignore')
from datetime import datetime

startTime = datetime.now()

def unstack_df(df):
    df_piv = pd.pivot(df, values='group', index=['chr','start', 'end'], columns=['query']).reset_index()
    col_df = df_piv.columns.values[:].tolist()
    return pd.DataFrame(df_piv.iloc[:, :].values, columns=col_df)

def file_to_list(file):
    file_l = pd.read_csv(file)['genotype'].tolist()
    return file_l

def select_genotypes(df, sample_list):
    filter_df = df.set_index(['chr', 'start', 'end']).filter(items=sample_list, axis=1)
    return filter_df.reset_index()

def transform_hap(df):
    colum_list=df.columns[3:].to_list()
    end = len(df)
    i=0
    dfs = []
    while i < end:
        df_row = df.iloc[i:i+1, :]
        duplicated_df = pd.concat([df_row]*int(df_row.iloc[:, 3:].max(axis=1)+1), ignore_index=True)
        duplicated_df.insert(3, 'haplotype', range(0, 1 + len(duplicated_df)-1), allow_duplicates=True)
        header = duplicated_df.iloc[:, :4]
        for sample in colum_list:
                header[sample] = np.where(duplicated_df[sample] == duplicated_df['haplotype'], 1, 0)
        dfs.append(header)
        i+=1
    header_concat = pd.concat(dfs, join="inner")
    header_concat.insert(3, 'position', header_concat['end'] + header_concat['haplotype'], allow_duplicates=True)
    header_concat['haplotype'] = header_concat['chr'] +'_'+ header_concat['position'].astype(str)
    return header_concat

def combine_chr(files, reference, chromosome, genotypes=None):
    # modify this part if multiple chromosomes are needed
    f_df = pd.read_csv(files, delimiter='\t')
    f_df_ref = f_df[(f_df['reference'] == reference) & (f_df['chromosome'].str.contains(chromosome))]
    dfs=[]
    for index, row in f_df_ref.iterrows():
        in_db = pd.read_csv(row['path']+row['file'], delimiter='\t')
        in_db = in_db[(in_db['group']>=0)] # TODO, find afaster way to filter negatives
        AP_df = unstack_df(in_db)
        if genotypes is not None:
            AP_filtered = select_genotypes(AP_df, genotypes)
            haplotypes_gwas = transform_hap(AP_filtered)
        else:
            haplotypes_gwas = transform_hap(AP_df)
        dfs.append(haplotypes_gwas)
    haplotypes_gwas_cat = pd.concat(dfs, axis=0)
    return haplotypes_gwas_cat

def get_plink_format(df, assembly):
    chromosome = df['chr'].str.split('_', expand=True)[0].values
    # num_lst = list(np.arange(0, len(df['position']),1))
    # pref = 'Haplo'
    # rs = [pref + str(sub) for sub in num_lst]
    rs = df['haplotype'].values
    dic = {'rs':rs,
           'alleles':'A/G',
           'chrom':chromosome,
           'strand':'+',
           'assembly':assembly,
           'center':'NA',
           'protLSID':'NA',
           'assayLSID':'NA',
           'panel':'NA',
           'Qcode':'NA'
          }
    heder = pd.DataFrame(dic)
    heder.insert(3, 'pos', df['position'].values)
    df_m = pd.concat([heder, df.iloc[:, 5:].reset_index(drop=True)], axis=1)
    tmp_df = df_m.set_index(list(df_m.columns[:11].values))
    tmp_df = tmp_df.where(tmp_df == 0, 'A')
    out_df = tmp_df.where(tmp_df == 'A', 'G').reset_index()
    return out_df

def create_dic(pfx_file):
    with open(pfx_file) as f:
     rows = (line.rstrip().split('\t') for line in f)
     dic = {row[0]:row[1] for row in rows}
     return dic

def map_substring(s, dict_map):
    for key in dict_map.keys():
        if key in s:
            return dict_map[key]
    return s

def rename_chr(df, dict_map):
    ref_col = df['chrom'].apply(lambda x: map_substring(x, dict_map))
    df['chrom'] = ref_col
    # df.insert(3, 'chrom', ref_col)
    return df

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--files')
    parser.add_argument('-r', '--reference')
    parser.add_argument('-c', '--chromosome')
    parser.add_argument('-l', '--lines', default=None)
    parser.add_argument('-s', '--plink', default=None)
    parser.add_argument('-p', '--pfx_file', default=None)
    parser.add_argument('-o', '--output')
    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()

    if args.lines is not None:
        genotypes_list = file_to_list(args.lines)
    else:
        genotypes_list = args.lines
    
    out_df = combine_chr(args.files, args.reference, args.chromosome, genotypes_list)
    if args.plink is not None:
        out_df = get_plink_format(out_df, args.reference)
        names_dic = create_dic(args.pfx_file)
        out_df = rename_chr(out_df, names_dic)

    out_df.to_csv(args.output, sep='\t', index=False)
    print(datetime.now() - startTime)
if __name__ == '__main__':
    main()