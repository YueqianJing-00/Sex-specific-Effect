#import essential packages
import re
from functools import reduce
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm, rankdata
import numpy as np
import os
import sys
import pickle
import statsmodels.api as sm
from pyplink import PyPlink
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm

traits = sorted(pd.read_csv('../file/traits.txt',header=None).iloc[:,0])
replication_traits = sorted(pd.read_csv('../file/replication_traits.csv',sep='\t')['traits'].to_list())
traitgroup_dict = pickle.load(open('../file/traitgroup_dict.pkl', 'rb'))
traitgroups = pickle.load(open('../file/traitgroups.pkl', 'rb'))
path_to_rep = '../rep/'
path_to_trait = '../Primary_Summary_Statistics/'


def test_SSE_first(df):
    if len(df) == 0:
        return df
    df = df.iloc[:,3:]
    
    primary_result = pd.read_csv('{0}{1}/{1}.clumped_all.tsv'.format(path_to_trait, trait), sep = '\t')
    primary_result = primary_result[primary_result['SNP'].isin(df['SNP'])]
    primary_result['Direction'] = np.where(primary_result['beta_female'] > primary_result['beta_male'], 1, -1)
    primary_result = df[['SNP']].merge(primary_result, on='SNP')
    df = df.merge(primary_result[['SNP', 'locus','Direction']], on='SNP', how='left')
    df['z_diff'] = (df['Beta_female'] - df['Beta_male']) / (df['SE_female']**2 + df['SE_male']**2)**.5
    df['p_diff'] = 1 - norm.cdf(df['z_diff'] * df['Direction'])   
    df['p_diff_adj'] = multipletests(df['p_diff'], alpha = 0.3, method = 'fdr_bh')[1]
    df = df[['traits','locus','SNP','Beta_male','SE_male','P_male','Beta_female','SE_female','P_female','z_diff','p_diff','p_diff_adj']]
    df = df[df['p_diff_adj'] < 0.05]
    df.to_csv('{0}{1}/{1}.replication_SSE_first_class.tsv'.format(path_to_rep, trait), sep = '\t',index=False)
    
    return df

def test_SSE_second(df):
    
    if len(df) == 0:
        return df
    df = df.iloc[:,3:]
        
    primary_result = pd.read_csv('{0}{1}/{1}.clumped_all.tsv'.format(path_to_trait, trait), sep = '\t')
    primary_result = primary_result[primary_result['SNP'].isin(df['SNP'])]
    primary_result['Direction'] = np.where(primary_result['beta_female'] > primary_result['beta_male'], 1, -1)
    primary_result = df[['SNP']].merge(primary_result, on='SNP')
    df = df.merge(primary_result[['SNP', 'locus','Direction']], on='SNP', how='left')
    df['z_diff'] = (df['Beta_female'] - df['Beta_male']) / (df['SE_female']**2 + df['SE_male']**2)**.5
    df['p_diff'] = 1 - norm.cdf(df['z_diff'] * df['Direction'])
    df['p_diff_adj'] = multipletests(df['p_diff'], alpha = 0.3, method = 'fdr_bh')[1]
    df = df[['traits','locus','SNP','Beta_male','SE_male','P_male','Beta_female','SE_female','P_female','z_diff','p_diff','p_diff_adj']]
    df = df[df['p_diff_adj'] < 0.05]
    df.to_csv('{0}{1}/{1}.replication_SSE_second_class.tsv'.format(path_to_rep, trait), sep = '\t',index=False)
    
    return df

i = int(sys.argv[1]) - 1
trait = replication_traits[i]

df = pd.read_csv('{0}{1}/{1}.replication_first_stage.tsv'.format(path_to_trait, trait), sep = '\t')
df_SSE_first = df[(df['P_male'] < 5e-8) | (df['P_female'] < 5e-8)]
df_SSE_second = df[((df['P_male'] < 5e-5) | (df['P_female'] < 5e-5))]

test_SSE_first(df_SSE_first)
test_SSE_second(df_SSE_second)

print("done")
