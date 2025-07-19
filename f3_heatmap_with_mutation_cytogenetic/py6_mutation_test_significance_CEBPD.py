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
        p1_total = df[col][df.index.intersection(patient_c1)].shape[0]
        p1_sum = df[col][df.index.intersection(patient_c1)].sum()
        p2_total = df[col][df.index.intersection(patient_c2)].shape[0]
        p2_sum = df[col][df.index.intersection(patient_c2)].sum()
        s,p = stats.fisher_exact([[p1_sum,p1_total-p1_sum],[p2_sum,p2_total-p2_sum]])
        sig_summary.loc[col,'CEBPD_down_total'] = p1_total
        sig_summary.loc[col,'CEBPD_down_mutated'] = p1_sum
        sig_summary.loc[col,'CEBPD_NotChanged_total'] = p2_total
        sig_summary.loc[col,'CEBPD_NotChanged_mutated'] = p2_sum
        sig_summary.loc[col,'Odds Ratio'] = s
        sig_summary.loc[col,'pvalue'] = p    
        sig_summary['fdr'] = fdr_adj_p(sig_summary['pvalue'],'pvalue')
    return sig_summary


### main section
indir = 'f6_indir_mutation_fisher_test_CEBPD'
outdir='f6_mutation_fisher_test_CEBPD'
os.makedirs(outdir,exist_ok=True)

project_dir="/nv/vol190/zanglab/zw5j/wwork2018/AML_fran"
project_dir="/Volumes/zanglab/zw5j/wwork2018/AML_fran"


# ==== mutation info, genes with # mutation >10 across all PAML patients, not limited to 59 paired
mutation_df = pd.read_csv('{}/ff5_CNV_mutation/f2_mutation/PAML_gene_mutation_relapse_GT10.csv'.format(project_dir),index_col=0)
mutation_df = mutation_df.fillna(0)
selected_genes = ['DNMT3A','FLT3_ITD','NPM1','NRAS','BRAF',
                  'SETBP1','TET2','WT1','ASXL1','CEBPA',
                  'PTPN11','FLT3','TP53','GATA2','KMT2D']
# selected_genes = mutation_df.sum().sort_values(ascending=False)[:15].index
mutation_df = mutation_df[selected_genes]


clinical_dir = '{}/f0_fran_data_new/data_processed/'.format(project_dir)
id_info = pd.read_csv(clinical_dir+os.sep+'Combined_clinical.csv',index_col=0)

# ==== cytogenetics 
# cytogenetic_df = pd.read_excel('{}/ff5_CNV_mutation/data/PAML_cytogenetics_category_col_F-S.xlsx'.format(project_dir),index_col=0)
# selected_cols=['NCG','Other','Complex']
# selected_cols=['NCG','Complex','Inv(16)','8+','t(8;21)','del5(q)','7-','17-/17p abn','11q23','5-','t(9;22)','t(6;9)','Inv(3) or t(3;3)']
# cytogenetic_df = cytogenetic_df.fillna(0)
# cytogenetic_df = cytogenetic_df[selected_cols]


# clinical info
# clinical_dir = '{}/f0_fran_data_new/data_processed/'.format(project_dir)

writer = pd.ExcelWriter(outdir+os.sep+'mutation_CEBPD_down_vs_Not_changed.xlsx')

basenames = ['PAML','SG','Combined']
# fisher exact test
for basename in basenames[:]:
    
    patient_df = pd.read_excel('{}/DEG_CEBPD_relAML.xlsx'.format(indir),sheet_name=basename,index_col=0)
    
    # patient_down = patient_df.dropna().index
    patient_down = patient_df[~patient_df.isna().CEBPD].index
    patient_nan = patient_df[patient_df.isna().CEBPD].index
    
    patient_c1 = id_info.loc[patient_down]['subject_ID']
    patient_c2 = id_info.loc[patient_nan]['subject_ID']
    
    sig_summary = return_sig(mutation_df,patient_c1,patient_c2)
    sig_summary = sig_summary.sort_values(by='pvalue',ascending=True)
    sig_summary.to_excel(writer,basename)

    # sig_summary = return_sig(cytogenetic_df,patient_c1,patient_c2)
    # sig_summary.to_csv(outdir+os.sep+'{}_cytogenetic_sig.csv'.format(basename))
writer.close()





