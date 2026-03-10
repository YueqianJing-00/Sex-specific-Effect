import sys
import pickle
import pandas as pd
from scipy.stats import norm

trait = sys.argv[1]

path = './Primary_Summary_Statistics/'

variant_dict = pickle.load(open('./variants.pkl', 'rb'))

all_df = pd.read_csv('{0}{1}/{1}.gwas.imputed_v3.both_sexes.tsv.bgz'.format(path, trait), sep = '\t', compression = 'gzip')
male_df = pd.read_csv('{0}{1}/{1}.gwas.imputed_v3.male.tsv.bgz'.format(path, trait), sep = '\t', compression = 'gzip')
female_df = pd.read_csv('{0}{1}/{1}.gwas.imputed_v3.female.tsv.bgz'.format(path, trait), sep = '\t', compression = 'gzip')

all_df['SNP'] = all_df['variant'].map(variant_dict)
male_df['SNP'] = male_df['variant'].map(variant_dict)
female_df['SNP'] = female_df['variant'].map(variant_dict)

all_df = all_df.rename(columns = {'pval': 'P', 'n_complete_samples': 'N', 'minor_allele': 'A1'})
male_df = male_df.rename(columns = {'pval': 'P', 'n_complete_samples': 'N', 'minor_allele': 'A1'})
female_df = female_df.rename(columns = {'pval': 'P', 'n_complete_samples': 'N', 'minor_allele': 'A1'})

all_df['A2'] = all_df['variant'].apply(lambda x: x.split(':')[2])
male_df['A2'] = male_df['variant'].apply(lambda x: x.split(':')[2])
female_df['A2'] = female_df['variant'].apply(lambda x: x.split(':')[2])

all_df.to_csv('{0}{1}/{1}.sumstats.both_sexes.tsv'.format(path, trait), sep = '\t', index = False, na_rep = 'NA')
male_df.to_csv('{0}{1}/{1}.sumstats.male.tsv'.format(path, trait), sep = '\t', index = False, na_rep = 'NA')
female_df.to_csv('{0}{1}/{1}.sumstats.female.tsv'.format(path, trait), sep = '\t', index = False, na_rep = 'NA')
