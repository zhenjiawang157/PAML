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
    # df = df[df<0]
    # df = np.abs(df)
    df = df.notnull()
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
outdir='f6b_mutation_changes_test'
os.makedirs(outdir,exist_ok=True)

project_dir="/nv/vol190/zanglab/zw5j/wwork2018/AML_fran"
project_dir="/Volumes/zanglab/zw5j/wwork2018/AML_fran"


# cohort I
basename = 'PAML_pthre5_rk3_ck2_patient'
mutation = pd.read_csv('data_mutation_CNV/mutation_DvsR_cohortI_matrix.csv',index_col=0)   
patient_c1 = [i.strip() for i in open(indir+'/{}1.txt'.format(basename)).readlines()]
patient_c2 = [i.strip() for i in open(indir+'/{}2.txt'.format(basename)).readlines()]
    
sig_summary = return_sig(mutation,patient_c1,patient_c2)
sig_summary = sig_summary.sort_values(by = 'pvalue')
# sig_summary.to_csv(outdir+os.sep+'{}_mutation_changes_cohortI.csv'.format(basename))


cnv = pd.read_csv('data_mutation_CNV/CNV_DvsR_cohortI_matrix.csv',index_col=0)   
cnv = cnv.loc[set(patient_c1).union(patient_c2)]
cnv = cnv.loc[:,cnv.notnull().sum(axis=0)>0]
sig_summary = return_sig(cnv,patient_c1,patient_c2)
sig_summary = sig_summary.sort_values(by = 'pvalue')
sig_summary.to_csv(outdir+os.sep+'{}_CNV_cohortI.csv'.format(basename))


# cohort II
basename = 'SG_pthre5_rk3_ck2_patient'
mutation = pd.read_csv('data_mutation_CNV/mutation_DvsR_cohortII_matrix.csv',index_col=0)   
patient_c1 = [i.strip() for i in open(indir+'/{}1.txt'.format(basename)).readlines()]
patient_c2 = [i.strip() for i in open(indir+'/{}2.txt'.format(basename)).readlines()]
    
sig_summary = return_sig(mutation,patient_c1,patient_c2)
sig_summary = sig_summary.sort_values(by = 'pvalue')
# sig_summary.to_csv(outdir+os.sep+'{}_mutation_changes_cohortII.csv'.format(basename))




