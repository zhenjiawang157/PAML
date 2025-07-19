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

def return_match_name(col,cluster_id):
    if col=='epigenetic_cluster':
        if cluster_id==1:
            return 'Diagnosis dominant'
        if cluster_id==2:
            return 'Mixed'
        if cluster_id==3:
            return 'Relapse dominant'
    
    if col=='genetic_cluster_new':
        if cluster_id==1:
            return 'Clonal changes'
        if cluster_id==2:
            return 'Stable'
        if cluster_id==3:
            return 'Subclonal changes'


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

    
    
def return_sig(indf,patient_c1,patient_c2,col):
    # df = df[df<0]
    # df = np.abs(df.col==cluster_id)
    sig_summary=pd.DataFrame()
    for cluster_id in [1,2,3]:
        df = indf[indf[col]==cluster_id]
        p1_total = df[col][patient_c1].shape[0]
        p1_sum = df[col][patient_c1].notnull().sum()
        p2_total = df[col][patient_c2].shape[0]
        p2_sum = df[col][patient_c2].notnull().sum()
        s,p = stats.fisher_exact([[p1_sum,p1_total-p1_sum],[p2_sum,p2_total-p2_sum]])
        sig_summary.loc[return_match_name(col,cluster_id),'p1_total'] = p1_total
        sig_summary.loc[return_match_name(col,cluster_id),'p1_sum'] = p1_sum
        sig_summary.loc[return_match_name(col,cluster_id),'p2_total'] = p2_total
        sig_summary.loc[return_match_name(col,cluster_id),'p2_sum'] = p2_sum
        sig_summary.loc[return_match_name(col,cluster_id),'stats_score'] = s
        sig_summary.loc[return_match_name(col,cluster_id),'pvalue'] = p    
        sig_summary['fdr'] = fdr_adj_p(sig_summary['pvalue'],'pvalue')
    return sig_summary


### main section
indir = 'f6_indir_patient_clusters_cp'
outdir='f6c_evolution_test'
os.makedirs(outdir,exist_ok=True)

project_dir="/nv/vol190/zanglab/zw5j/wwork2018/AML_fran"
project_dir="/Volumes/zanglab/zw5j/wwork2018/AML_fran"


# cohort I
basename = 'PAML_pthre5_rk3_ck2_patient'
evolution = pd.read_csv('data_mutation_CNV/evolution_corhotI.csv',index_col=0)   
patient_c1 = [i.strip() for i in open(indir+'/{}1.txt'.format(basename)).readlines()]
patient_c2 = [i.strip() for i in open(indir+'/{}2.txt'.format(basename)).readlines()]

for cluster_col in ['epigenetic_cluster','genetic_cluster_new']:
    sig_summary = return_sig(evolution,patient_c1,patient_c2,cluster_col)
    sig_summary = sig_summary.sort_values(by = 'pvalue')
    sig_summary.to_csv(outdir+os.sep+'{}_{}_cohortI.csv'.format(basename,cluster_col))


# cohort II
basename = 'SG_pthre5_rk3_ck2_patient'
evolution = pd.read_csv('data_mutation_CNV/evolution_corhotII.csv',index_col=0)   
patient_c1 = [i.strip() for i in open(indir+'/{}1.txt'.format(basename)).readlines()]
patient_c2 = [i.strip() for i in open(indir+'/{}2.txt'.format(basename)).readlines()]
    
for cluster_col in ['epigenetic_cluster']:
    sig_summary = return_sig(evolution,patient_c1,patient_c2,cluster_col)
    sig_summary = sig_summary.sort_values(by = 'pvalue')
    sig_summary.to_csv(outdir+os.sep+'{}_{}_cohortII.csv'.format(basename,cluster_col))




