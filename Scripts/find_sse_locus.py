import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import sys
import pickle
from tqdm import tqdm
from statsmodels.stats.multitest import multipletests

traitgroup_num = int(sys.argv[1])

def classify_effect_direction(group):
    total = len(group)
    consistent = sum((group['beta_male'] > 0) == (group['beta_female'] > 0))
    if consistent/total >= 0.8:
        return 'CED'
    elif consistent/total <= 0.2:
        return 'OED'
    else:
        return 'AED'

traits = pd.read_csv('../files/traits.txt')['Trait'].to_list()
traitgroup_dict = pickle.load(open('../files/traitgroup_dict.pkl', 'rb'))
traitgroups = pickle.load(open('../files/traitgroups.pkl', 'rb'))

for i in range(len(traitgroups[traitgroup_num - 1])):
    trait = traitgroups[traitgroup_num - 1][i]
    SSE = pd.read_csv('../Primary_Summary_Statistics/{0}/{0}.clumped_SSE.tsv'.format(trait), sep = '\t')
    SSE = SSE[SSE['p_diff_adj'] < 0.05]
    non_SSE = pd.read_csv('../Primary_Summary_Statistics/{0}/{0}.clumped_all.tsv'.format(trait), sep = '\t')

    SSE_representative = SSE.loc[SSE.groupby('locus')['p_min'].idxmin()]
    SSE_representative.to_csv('../Primary_Summary_Statistics/{0}/{0}.clumped_SSE_representative.tsv'.format(trait), sep = '\t', index = False)
    
    #create a ED_candiadate dataframe that get rid of all the locus in SSE
    ED_candidate = non_SSE[~non_SSE['locus'].isin(SSE['locus'])]
    ED_candidate = ED_candidate[((ED_candidate['P_male'] < 5e-8) & (ED_candidate['P_female'] < 0.05))
              | ((ED_candidate['P_female'] < 5e-8) & (ED_candidate['P_male'] < 0.05))]
    
    if ED_candidate.shape[0] > 0:
        Direction = ED_candidate.groupby('locus').apply(classify_effect_direction)
        ED_summary = ED_candidate.loc[ED_candidate.groupby('locus')['p_min'].idxmin()]
        ED_summary['Direction'] = Direction.values
        ED_summary.to_csv('../Primary_Summary_Statistics/{0}/{0}.clumped_ED_summary.tsv'.format(trait), sep = '\t', index = False)
