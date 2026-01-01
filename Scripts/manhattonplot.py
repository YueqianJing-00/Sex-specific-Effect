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


def convert_position_to_ind_P(snp,sex):
    if sex == 'female':
        df_tmp = df_female[df_female['SNP']==snp]
    else:
        df_tmp = df_male[df_male['SNP']==snp]
    
    return df_tmp['ind'].tolist()[0]

traitgroup_num = int(sys.argv[1])
traits = sorted(pd.read_csv('../files/traits.txt',header=None).iloc[:,0])
replication_traits = sorted(pd.read_csv('../files/replication_traits.csv',sep='\t')['traits'].to_list())
traitgroup_dict = pickle.load(open('../files/traitgroup_dict.pkl', 'rb'))
traitgroups = pickle.load(open('../files/traitgroups.pkl', 'rb'))
path_to_trait = '../Primary_Summary_Statistics/'

variant = pd.read_csv("../files/snp_ref_alt_pos_chr.txt",sep='\t',dtype={"chr": str})

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


locus_list = rep_groups[traitgroup_num]
for snp in locus_list:
    snp_index.append(int(variant[variant.SNP == snp].index[0]))
snp_index.sort()
locus_list = variant.iloc[snp_index,:].SNP.tolist()

pval_list_male = []
pval_list_female = []
snp_list_male = []
snp_list_female = []

for locus in tqdm(locus_list):
    tmp_male = []
    tmp_female = []
    tmp_snp_male = []
    tmp_snp_female = []
    for i in range(len(traitgroups[traitgroup_num])):
        trait = traitgroups[traitgroup_num][i]
        try:
            primary_sig_sse = pd.read_csv("{0}{1}/{1}.clumped_all.tsv".format(path_to_trait,trait),sep='\t')
            pval_locus = primary_sig_sse[primary_sig_sse.locus == locus]
            p_male = pval_locus.P_male.tolist()
            p_female = pval_locus.P_female.tolist()
            snp = pval_locus.SNP.tolist()
            if len(pval_locus) > 0:
                tmp_male.append(min(p_male))
                tmp_female.append(min(p_female))
                tmp_snp_male.append(snp[np.argmin(np.array(p_male))])
                tmp_snp_female.append(snp[np.argmin(np.array(p_female))])
        except:
            continue
    
    pval_list_male.append(min(tmp_male))
    pval_list_female.append(min(tmp_female))
    snp_list_male.append(tmp_snp_male[np.argmin(np.array(tmp_male))])
    snp_list_female.append(tmp_snp_female[np.argmin(np.array(tmp_female))])


fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

x_labels = []
x_labels_pos = []
y_max = 0

texts_female = []
texts_male = []
    
for i in tqdm(range(len(traitgroups[traitgroup_num]))):
    trait = traitgroups[traitgroup_num][i]
    df = pd.read_csv("{0}{1}/{1}.sumstats.tsv".format(path_to_trait,trait),sep='\t')
    
    df_male = pd.DataFrame(columns = ['chr','BP','SNP','P'])
    df_male['chr'] = variant.chr
    df_male['BP']= variant.pos
    df_male['SNP']= variant.SNP
    df_male.P = df.P_male
    df_male.P = df_male.P.apply(lambda x: np.nan if x <= 0 else x)
    #pval_matrix_male = pd.concat([pval_matrix_male,pd.DataFrame(df_male.iloc[snp_index,3])],axis=1)
    #df_male.replace([np.inf, -np.inf], np.nan, inplace=True)
    df_male.dropna(subset=['P'],inplace=True)
    df_male['P'] = np.log10(df_male.P)
    df_male['ind'] = range(len(df_male))
    df_grouped_male = df_male.groupby(('chr'))
    
    df_female = pd.DataFrame(columns = ['chr','BP','SNP','P'])
    df_female['chr'] = variant.chr
    df_female['BP']= variant.pos
    df_female['SNP']= variant.SNP
    df_female.P = df.P_female
    df_female.P = df_female.P.apply(lambda x: np.nan if x <= 0 else x)
    #pval_matrix_female = pd.concat([pval_matrix_female,pd.DataFrame(df_female.iloc[snp_index,3])],axis=1)
    #df_male.replace([np.inf, -np.inf], np.nan, inplace=True)
    df_female.dropna(subset=['P'],inplace=True)
    df_female['P'] = -np.log10(df_female.P)
    df_female['ind'] = range(len(df_female))
    df_grouped_female = df_female.groupby(('chr'))
    
    y_max_tmp = max(-min(df_male.P),max(df_female.P))
    y_max_tmp = math.ceil(y_max_tmp / 10) * 10
    y_max = max(y_max, y_max_tmp)

# Create a figure with two subplots sharing the x-axis

    # The first plot (this will be on the top)
    for num, (name, group) in enumerate(df_grouped_female):
        if num == 0:
            ax1.scatter(x=group['ind'], y=group['P'], color=colors[len(traitgroups[traitgroup_num])-i-1], s=0.8,label=trait)
        else:
            ax1.scatter(x=group['ind'], y=group['P'], color=colors[len(traitgroups[traitgroup_num])-i-1], s=0.8)
        if i == 0:
            if num > 0:
                ax1.axvline(x=group['ind'].iloc[0], color='grey', linestyle='-', linewidth=0.75)

    # The second plot (this will be on the bottom, inverted)
    for num, (name, group) in enumerate(df_grouped_male):
        if num == 0:
            ax2.scatter(x=group['ind'], y=group['P'], color=colors[len(traitgroups[traitgroup_num])-i-1], s=0.8,label=trait)
        else:
            ax2.scatter(x=group['ind'], y=group['P'], color=colors[len(traitgroups[traitgroup_num])-i-1], s=0.8)
        if i == 0:
            x_labels.append(name)
            x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))
            if num > 0:
                ax2.axvline(x=group['ind'].iloc[0], color='grey', linestyle='-', linewidth=0.75)

if len(pval_list_male) > 0:
    y_max = min(y_max, 6 * np.max(-np.log10(pval_list_male)), 6 * np.max(-np.log10(pval_list_female)))


# Set the spines to zero for both subplots
ax1.spines['left'].set_position('zero')
#ax1.spines['right'].set_color('none')
#ax1.spines['top'].set_color('none')
ax1.set_xlim([0, len(df)])
ax1.set_ylim(bottom=0,top=y_max)  # Set a uniform ylim for both, adjust 'top' as needed

ax2.spines['left'].set_position('zero')
#ax2.spines['right'].set_color('none')
#ax2.spines['top'].set_color('none')
ax2.set_xlim([0, len(df)])
ax2.set_ylim(top=0,bottom=-y_max)  # Set a uniform ylim for both, adjust 'top' as needed

# Rotate the x-axis labels to 45 degrees on the bottom subplot
ax2.set_xticks(x_labels_pos)
ax2.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=7.5)

# Set common x and y labels
ax2.set_xlabel('Chromosome')
ax1.set_ylabel('Female -log(p) value')
ax2.set_ylabel('Male -log(p) value')

ax1.axhline(y= -np.log10(5e-8), color='grey', linestyle='-', linewidth=1)
ax2.axhline(y= np.log10(5e-8), color='grey', linestyle='-', linewidth=1)

# Adjust the layout
plt.tight_layout()
plt.subplots_adjust(right=0.9)

# Adjust the spacing between the subplots
plt.subplots_adjust(hspace=0)  # hspace=0 makes no space between subplots

# Remove the x-axis tick marks from the top subplot
ax1.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

flag = 0
snp_num = len(rep_groups[traitgroup_num])

#pval_list = np.min(pval_matrix,axis=1).tolist()
pval_list_male = -np.log10(pval_list_male)
pval_list_female = -np.log10(pval_list_female)

for i in tqdm(range(snp_num)) :
    locus = locus_list[i]
    info = nearby_gene_groups[locus]
    # Convert chromosomal position to 'ind' if necessary, depending on how your df is structured
    
    text = ', '.join(info['gene'])
    
    if pval_list_female[i] > pval_list_male[i]:
        snp = snp_list_female[i]
        ind = convert_position_to_ind_P(snp,'female')
        Pval = pval_list_female[i]
        #texts_female.append(ax1.text(ind,  Pval, text, ha='right', va='top'))
    # Annotate on the top subplot, adjust the y-value as needed
        ax1.annotate(text, xy=(ind, Pval), xytext=(int((flag+1) * len(df_female)/(snp_num+1)) , (y_max-1) * ( 1 - 0.05 * (flag % 6))  ),
                 arrowprops=dict(facecolor='black', arrowstyle='->', lw=0.5),
                 horizontalalignment='right', verticalalignment='top')
    else:
        snp = snp_list_male[i]
        ind = convert_position_to_ind_P(snp,'male')
        Pval = pval_list_male[i]
        #texts_male.append(ax2.text(ind, - Pval, text, ha='right', va='bottom'))
    # Annotate on the bottom subplot, adjust the y-value as needed
        ax2.annotate(text, xy=(ind, -Pval), xytext=(int((flag+1) * len(df_female)/(snp_num+1)), -(y_max-1) * ( 1 - 0.03 * (flag % 6))  ),
                 arrowprops=dict(facecolor='black', arrowstyle='->', lw=0.5),
                 horizontalalignment='right', verticalalignment='bottom')
    
    flag += 1

#adjust_text(texts_female, ax=ax1, arrowprops=dict(facecolor='black', arrowstyle='->', lw=1.5))
#adjust_text(texts_male, ax=ax2, arrowprops=dict(facecolor='black', arrowstyle='->', lw=1.5))    
handle, label = ax1.get_legend_handles_labels()
fig.legend(handle, label, loc='upper right', bbox_to_anchor=(1, 1))

print(traitgroup_num)
fig.savefig('../result/manhattan_plot/traitgroup_{0}.png'.format(traitgroup_num), dpi=500, bbox_inches='tight')
#plt.show()
