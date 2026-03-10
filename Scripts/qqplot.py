import pandas as pd
import numpy as np
import pickle
import sys
import os
from pyplink import PyPlink
from tqdm import tqdm
from math import floor,log10
from decimal import Decimal, getcontext, ROUND_HALF_EVEN
from scipy.stats import uniform
from scipy.stats import randint
import matplotlib.pyplot as plt
import math
import pickle
from matplotlib.lines import Line2D
import re
 
traitgroup_num = int(sys.argv[1])
print(traitgroup_num)

colors = [
    "#1f77b4",  # Muted blue
    "#ff7f0e",  # Safety orange
    "#2ca02c",  # Cooked asparagus green
    "#d62728",  # Brick red
    "#9467bd",  # Muted purple
    "#8c564b",  # Chestnut brown
    "#e377c2",  # Raspberry yogurt pink
    "#7f7f7f",  # Middle gray
    "#bcbd22",  # Curry yellow-green
    "#17becf",  # Blue-teal
    "#aec7e8",  # Soft blue
    "#ffbb78",  # Melon orange
]




traits = sorted(pd.read_csv('../files/traits.txt',header=None).iloc[:,0])
replication_traits = sorted(pd.read_csv('../files/replication_traits.csv',sep='\t')['traits'].to_list())
traitgroup_dict = pickle.load(open('../files/traitgroup_dict.pkl', 'rb'))
traitgroups = pickle.load(open('../files/traitgroups.pkl', 'rb'))
path_to_trait = '../Primary_Summary_Statistics/'

# Assuming 'traits' is a list of trait identifiers
# Initialize an empty DataFrame outside the loop
results_df = pd.DataFrame()

for trait in tqdm(traits):
    try:
        with open('{0}{1}/{1}_ldsc_result_male_female.log'.format(path_to_trait,trait), 'r') as file:
            log_content = file.read()

        # Define the regular expression patterns for lambda GC and intercept
        lambda_gc_pattern = r'Lambda GC: ([\d\.]+)'
        intercept_pattern = r'Intercept: ([\d\.]+) \(([\d\.]+)\)'

        # Find all matches for lambda GC and intercept
        lambda_gc_matches = re.findall(lambda_gc_pattern, log_content)
        intercept_matches = re.findall(intercept_pattern, log_content)

        # Extract the values for lambda GC and intercept for both traits
        lambda_gc_values = [float(match) for match in lambda_gc_matches]
        intercept_values = [(float(match[0]), float(match[1])) for match in intercept_matches]

        # Create a temporary DataFrame for the current trait
        df_tmp = pd.DataFrame({
            'trait': [trait],
            'lambda_male': [lambda_gc_values[0]],
            'lambda_female': [lambda_gc_values[1]],
            'intercept_male': [intercept_values[0][0]],
            'intercept_female': [intercept_values[1][0]],
            'intercept_male_se': [intercept_values[0][1]],
            'intercept_female_se': [intercept_values[1][1]]
        }, index=[0])                   

        # Append the temporary DataFrame to the main results DataFrame
        results_df = pd.concat([results_df, df_tmp], ignore_index=True)
    except:
        print(trait)

plt.figure(figsize=(12,9))

legend_elements = []
male_markers = []
female_markers = []
labels = []
legend_handles = []

for i in tqdm(range(len(traitgroups[traitgroup_num]))):
    trait = traitgroups[traitgroup_num][i]
    df = pd.read_csv("{0}{1}/{1}.sumstats.tsv".format(path_to_trait,trait),sep='\t')
    
    P_male = np.array(df.P_male.tolist())
    P_female = np.array(df.P_female.tolist())
    
    P_male.sort()
    P_female.sort()
    
    P_male = P_male[~np.isnan(P_male)]
    P_female = P_female[~np.isnan(P_female)]
    
    expected_P_male = [(i - 0.5) / len(P_male) for i in range(1, len(P_male) + 1)]
    expected_P_female = [(i - 0.5) / len(P_female) for i in range(1, len(P_female) + 1)]


    # Convert to -log10 scale
    expected_P_male = -np.log10(expected_P_male)
    P_male = -np.log10(P_male)

    expected_P_female = -np.log10(expected_P_female)
    P_female = -np.log10(P_female)
    
    color = colors[len(traitgroups[traitgroup_num])-i-1]
    #labels.append(trait)
    legend_elements.append(plt.scatter([], [], s=50, facecolors='none', edgecolors=color, marker='^'))
    legend_elements.append(plt.scatter([], [], s=50, facecolors='none', edgecolors=color, marker='o'))
    #legend_elements.append(male_markers)
    #legend_elements.append(female_markers)
    try:
        intercept_male = results_df.loc[results_df['trait'] == trait, 'intercept_male'].values[0]
        lambda_gc_male = results_df.loc[results_df['trait'] == trait, 'lambda_male'].values[0]
        intercept_female = results_df.loc[results_df['trait'] == trait, 'intercept_female'].values[0]
        lambda_gc_female = results_df.loc[results_df['trait'] == trait, 'lambda_female'].values[0]
        
        print(intercept_male)
        print(intercept_female)

        # Create the Q-Q plot
        plt.scatter(expected_P_male, P_male, s = 3.5, facecolors='none', edgecolors = colors[len(traitgroups[traitgroup_num])-i-1], marker='^')
        plt.scatter(expected_P_female, P_female, s = 3.5 , facecolors='none', edgecolors = colors[len(traitgroups[traitgroup_num])-i-1], marker='o')

        male_line = Line2D([0], [0], marker='^', color='none', markerfacecolor='none',
                           markeredgecolor=color, markersize=3, label=f'{trait} Male Intercept: {intercept_male:.3f}, λ_GC: {lambda_gc_male:.3f}')
        female_line = Line2D([0], [0], marker='o', color='none', markerfacecolor='none',
                             markeredgecolor=color, markersize=3, label=f'{trait} Female Intercept: {intercept_female:.3f}, λ_GC: {lambda_gc_female:.3f}')

        # Add the handles to the list
        legend_handles.extend([male_line, female_line])
    
    except:
        plt.scatter(expected_P_male, P_male, s = 3.5, facecolors='none', edgecolors = colors[len(traitgroups[traitgroup_num])-i-1], marker='^')
        plt.scatter(expected_P_female, P_female, s = 3.5 , facecolors='none', edgecolors = colors[len(traitgroups[traitgroup_num])-i-1], marker='o')

        male_line = Line2D([0], [0], marker='^', color='none', markerfacecolor='none',
                           markeredgecolor=color, markersize=3, label=f'{trait} Male')
        female_line = Line2D([0], [0], marker='o', color='none', markerfacecolor='none',
                             markeredgecolor=color, markersize=3, label=f'{trait} Female')

        # Add the handles to the list
        legend_handles.extend([male_line, female_line])
    


# Add axis labels
plt.plot([0, max(expected_P_male)], [0, max(expected_P_male)], color='black', ls='-')  # Diagonal line
plt.xlabel('Expected -log10(p-value)')
plt.ylabel('Observed -log10(p-value)')

# Add a title
plt.title('Q-Q Plot of GWAS p-values')


# Construct the legend using the handles
plt.legend(handles=legend_handles, loc='upper left', markerscale=1, scatterpoints=1, fontsize=8)
plt.savefig('../result/qqplot/traitgroup_{0}.png'.format(traitgroup_num), dpi=500, bbox_inches='tight')
    
