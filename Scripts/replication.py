#import essential packages
import re
from functools import reduce
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm, rankdata
import os
import sys
import pickle
import statsmodels.api as sm
from pyplink import PyPlink
from tqdm import tqdm
# Function Definitions

i = int(sys.argv[1]) - 1

traits = sorted(pd.read_csv('../files/traits.txt',header=None).iloc[:,0])
replication_traits = sorted(pd.read_csv('../files/replication_traits.csv',sep='\t')['traits'].to_list())
traitgroup_dict = pickle.load(open('../files/traitgroup_dict.pkl', 'rb'))
traitgroups = pickle.load(open('../files/traitgroups.pkl', 'rb'))
path = '../files/'
path_to_trait = '../Primary_Summary_Statistics/'

#Read Covariate Data 
cov_both = pd.read_csv("../files/phen/cov_both.tsv",sep='\t')
cov_male = pd.read_csv("../files/phen/cov_male.tsv",sep='\t')
cov_female = pd.read_csv("../files/phen/cov_female.tsv",sep='\t')

def summary_model(model):
    coef = model.params['genotype']
    se = model.bse['genotype']
    p_value = model.pvalues['genotype']
    
    summary = []
    summary.append(coef)
    summary.append(se)
    summary.append(p_value)
    
    return summary

def linear_regression_gwas(trait):
    
    bed =  PyPlink('{0}{1}/{1}'.format(path_to_trait, trait))
    
    #Read Phenotypes Data
    phen_both = pd.read_csv('{0}{1}/phen_both.tsv'.format(path_to_trait, trait), sep = '\t',header = None)
    phen_female = pd.read_csv('{0}{1}/phen_female.tsv'.format(path_to_trait, trait), sep = '\t',header = None)
    phen_male = pd.read_csv('{0}{1}/phen_male.tsv'.format(path_to_trait, trait), sep = '\t',header = None)

    phen_both.columns = ['Individual ID','Phenotype']
    phen_female.columns = ['Individual ID','Phenotype']
    phen_male.columns = ['Individual ID','Phenotype']

    df = pd.DataFrame()

    for snp_name, genotypes in bed:

        genotypes = pd.DataFrame(genotypes)
        genotypes.index = genotypes.index 
        genotypes.columns = ['genotype']

        cov_1 = cov_both.iloc[:,2:]
        predict_both = pd.concat([genotypes,cov_1], axis=1)
        predict_both = sm.add_constant(predict_both)

        response_both = phen_both.iloc[:,1]
        response_both = pd.DataFrame(response_both)
        response_both.index = response_both.index 

        NaN_sample = response_both.isna().any(axis=1)
        response_both = response_both[-NaN_sample]
        predict_both = predict_both[-NaN_sample]

        male_sample = predict_both['sex'] == 1
        female_sample = predict_both['sex'] == 2

        predict_male = predict_both[male_sample]
        predict_male = predict_male.drop(columns = ['sex','sex_age'])
        response_male = response_both[male_sample]

        predict_female = predict_both[female_sample]
        predict_female = predict_female.drop(columns = ['sex','sex_age'])
        response_female = response_both[female_sample]

        model_both = sm.OLS(response_both, predict_both).fit()
        model_female = sm.OLS(response_female, predict_female).fit()
        model_male = sm.OLS(response_male, predict_male).fit()

        # Get summary statistics
        summary_both = summary_model(model_both)
        summary_female = summary_model(model_female)
        summary_male = summary_model(model_male)

        summary_all = summary_both + summary_female + summary_male
        summary_all = pd.DataFrame([summary_all], columns=[f'col{i}' for i in range(len(summary_all))])
        summary_all.columns = ['Beta_both','SE_both','P_both','Beta_female','SE_female','P_female','Beta_male','SE_male','P_male']
        summary_all['SNP'] = snp_name
        summary_all['traits'] = trait

        df = pd.concat([df, summary_all], axis=0)

    df.reset_index(drop=True, inplace=True)
    df.index = df.index + 1
    
    df.to_csv('{0}{1}/{1}.replication_first_stage.tsv'.format(path_to_trait, trait), sep = '\t',index=False)

    return 


trait = replication_traits[i]
linear_regression_gwas(trait)

print("done")
