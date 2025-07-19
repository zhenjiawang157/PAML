import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
#from GenomeData import *
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
# matplotlib.rcParams['font.size']=9
# import seaborn as sns
# sns.set(font_scale=1)
# sns.set_style("whitegrid", {'axes.grid' : False})
# from matplotlib import gridspec
#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')
# matplotlib.rcParams["font.sans-serif"] = ["Arial"]
# from  scipy.cluster.hierarchy import fcluster
from scipy import stats


def fdr_adj_p(pvalues,p_index):
    
    df = pvalues.to_frame()
    n,k = len(pvalues),len(pvalues)      
    minimum = 1    
    for index in df.sort_values(by=p_index,ascending=False).index:
        pvalue = df.loc[index,p_index]
        fdr = n*pvalue/k  
        minimum = min(minimum,fdr)
        df.loc[index,'fdr'] = minimum
        k=k-1     
    return df['fdr']


def return_sig(df,patient_c1,patient_c2):
    sig_summary=pd.DataFrame()
    for col in df.columns:
        p1_total = df[col][patient_c1].shape[0]
        p1_sum = df[col][patient_c1].sum()
        p2_total = df[col][patient_c2].shape[0]
        p2_sum = df[col][patient_c2].sum()
        s,p = stats.fisher_exact([[p1_sum,p1_total-p1_sum],[p2_sum,p2_total-p2_sum]])
        sig_summary.loc[col,'p1_total'] = p1_total
        sig_summary.loc[col,'p1_sum'] = p1_sum
        sig_summary.loc[col,'p2_total'] = p2_total
        sig_summary.loc[col,'p2_sum'] = p2_sum
        sig_summary.loc[col,'stats_score'] = s
        sig_summary.loc[col,'pvalue'] = p    
        sig_summary['fdr'] = fdr_adj_p(sig_summary['pvalue'],'pvalue')
    return sig_summary


### main section
indir = 'f6_indir_patient_clusters_cp'
outdir='f6_mutation_cytogenetic_fisher_test'
os.makedirs(outdir,exist_ok=True)

project_dir="/nv/vol190/zanglab/zw5j/wwork2018/AML_fran"
project_dir="/Volumes/zanglab/zw5j/wwork2018/AML_fran"


# ==== mutation info, genes with # mutation >10 across all PAML patients, not limited to 59 paired
mutation_df = pd.read_csv('{}/ff5_CNV_mutation/f2_mutation/PAML_gene_mutation_relapse_GT10.csv'.format(project_dir),index_col=0)
mutation_df = mutation_df.fillna(0)
selected_genes = ['DNMT3A','FLT3_ITD','NPM1','NRAS','BRAF',\
                  'SETBP1','TET2','WT1','ASXL1','CEBPA']
mutation_df = mutation_df[selected_genes]

# ==== cytogenetics 
cytogenetic_df = pd.read_excel('{}/ff5_CNV_mutation/data/PAML_cytogenetics_category_col_F-S.xlsx'.format(project_dir),index_col=0)
selected_cols=['NCG','Other','Complex']
selected_cols=['NCG','Complex','Inv(16)','8+','t(8;21)','del5(q)','7-','17-/17p abn','11q23','5-','t(9;22)','t(6;9)','Inv(3) or t(3;3)']
cytogenetic_df = cytogenetic_df.fillna(0)
cytogenetic_df = cytogenetic_df[selected_cols]


# clinical info
clinical_dir = '{}/f0_fran_data_new/data_processed/'.format(project_dir)

basenames = ['PAML_pthre5_rk3_ck2_euclidean_ward_cluster','PAML_pthre5_rk3_ck2_patient','SG_pthre5_rk3_ck2_patient']
# fisher exact test
for basename in basenames[:]:
    
    patient_c1 = [i.strip() for i in open(indir+'/{}1.txt'.format(basename)).readlines()]
    patient_c2 = [i.strip() for i in open(indir+'/{}2.txt'.format(basename)).readlines()]
        
    sig_summary = return_sig(mutation_df,patient_c1,patient_c2)
    sig_summary.to_csv(outdir+os.sep+'{}_mutation_sig.csv'.format(basename))

    sig_summary = return_sig(cytogenetic_df,patient_c1,patient_c2)
    sig_summary.to_csv(outdir+os.sep+'{}_cytogenetic_sig.csv'.format(basename))






