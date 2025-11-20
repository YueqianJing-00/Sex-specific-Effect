import numpy as np
import pandas as pd
from tqdm import tqdm
import os
import sys
import pickle
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests

def classify_effect_direction(group):
    total = len(group)
    consistent = sum((group['beta_male'] > 0) == (group['beta_female'] > 0))
    if consistent/total >= 0.8:
        return 'CED'
    elif consistent/total <= 0.2:
        return 'OED'
    else:
        return 'AED'

traitgroup_num = int(sys.argv[1])

traits = pd.read_csv('/gpfs/gibbs/pi/zhao/yj348/sse/traits.txt')['Trait'].to_list()
traitgroup_dict = pickle.load(open('/gpfs/gibbs/pi/zhao/yj348/sse/traitgroup_dict.pkl', 'rb'))
traitgroups = pickle.load(open('/gpfs/gibbs/pi/zhao/yj348/sse/traitgroups.pkl', 'rb'))
path = '/gpfs/gibbs/pi/zhao/yj348/sse/'
variant_dict = pickle.load(open('/gpfs/gibbs/pi/zhao/yj348/sse/variants.pkl', 'rb'))


df = pd.read_csv('/gpfs/gibbs/pi/zhao/yj348/sse/phen/Trait_Groups/{0}/{0}.clumped'.format(traitgroup_num), delim_whitespace=True,engine='python')
df['SNPs'] = df['SP2'].apply(lambda x: [y[:-3] for y in x.split(',')] if x != 'NONE' else np.nan)

clump_dict = df.set_index('SNP').to_dict()['SNPs']
clump_dict = {k: ([k] if type(v) == float else v + [k]) for k,v in clump_dict.items()}
#trait_clump_df_CED = pd.DataFrame(columns = ['trait', 'locus', 'SNP', 'P_male', 'beta_male', 'se_male', 'P_female', 'beta_female', 'se_female', 'p_diff'])
#trait_clump_df_OED = pd.DataFrame(columns = ['trait', 'locus', 'SNP', 'P_male', 'beta_male', 'se_male', 'P_female', 'beta_female', 'se_female', 'p_diff'])
#trait_clump_df_AED = pd.DataFrame(columns = ['trait', 'locus', 'SNP', 'P_male', 'beta_male', 'se_male', 'P_female', 'beta_female', 'se_female', 'p_diff'])


for i in range(len(traitgroups[traitgroup_num - 1])):
    trait_df = pd.read_csv('/gpfs/gibbs/pi/zhao/yj348/sse/Primary_Summary_Statistics/{0}/{0}.sumstats.tsv'.format(traitgroups[traitgroup_num-1][i]), sep = '\t')
    trait_clump_df_SSE = pd.DataFrame(columns = ['trait', 'locus', 'SNP', 'P_male', 'beta_male', 'se_male', 'P_female', 'beta_female', 'se_female', 'z_diff', 'p_diff', 'p_min'])
    trait_clump_df_no_SSE = pd.DataFrame(columns = ['trait', 'locus', 'SNP', 'P_male', 'beta_male', 'se_male', 'P_female', 'beta_female', 'se_female', 'z_diff', 'p_diff', 'p_min'])
    trait_clump_df = pd.DataFrame(columns = ['locus', 'SNP', 'P_male', 'beta_male', 'se_male', 'P_female', 'beta_female', 'se_female', 'z_diff', 'p_diff', 'p_min'])
    for locus, snps in tqdm(clump_dict.items()):
        locus_df = trait_df[trait_df['SNP'].isin(snps)].sort_values(by = 'p_min')  # sort by p_diff to get the smallest p_diff 
        locus_df = locus_df[['SNP', 'P_male', 'beta_male', 'se_male', 'P_female', 'beta_female', 'se_female', 'z_diff','p_diff', 'p_min']]
        locus_df['locus'] = locus
        trait_clump_df = pd.concat([trait_clump_df, locus_df])
    trait_clump_df_SSE = trait_clump_df[((trait_clump_df['P_male'] < 5e-8) & (trait_clump_df['P_female'] > 5e-8))
              | ((trait_clump_df['P_female'] < 5e-8) & (trait_clump_df['P_male'] > 5e-8))]
    trait_clump_df_no_SSE = trait_clump_df
    if len(trait_clump_df_SSE) != 0:
        trait_clump_df_SSE['p_diff_adj'] = multipletests(trait_clump_df_SSE['p_diff'], method = 'fdr_bh')[1]

    print("Done")
    print(traitgroups[traitgroup_num-1][i])

    #### Trait_clump_df_SSE contains all Sex-specific SNP grouped by locus
    trait_clump_df_SSE.to_csv('/gpfs/gibbs/pi/zhao/yj348/sse/Primary_Summary_Statistics/{0}/{0}.clumped_SSE.tsv'.format(traitgroups[traitgroup_num-1][i]), sep = '\t', index = False)
    #### Trait_clump_df_no_SSE contains all SNPs grouped by locus
    trait_clump_df_no_SSE.to_csv('/gpfs/gibbs/pi/zhao/yj348/sse/Primary_Summary_Statistics/{0}/{0}.clumped_all.tsv'.format(traitgroups[traitgroup_num-1][i]), sep = '\t', index = False)

print("Finish")
