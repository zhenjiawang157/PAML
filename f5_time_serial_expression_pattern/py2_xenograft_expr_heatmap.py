import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
from scipy import stats
#from GenomeData import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=12
import seaborn as sns
sns.set(font_scale=1)
sns.set_style("whitegrid", {'axes.grid' : False})
from matplotlib import gridspec
#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
from  scipy.cluster.hierarchy import fcluster


def quantile_normalization(dataframe):
    '''
    dataframe with samples in columns and probes across rows
    '''
    df = dataframe
    rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
    df_normalized=df.rank(method='min').stack().astype(int).map(rank_mean).unstack()
    return df_normalized

def log_quantile_normalization(df):
    '''
    dataframe with samples in columns and probes across rows
    '''
    df = np.log2(df+1)
    rank_mean = df.stack().groupby(df.rank(method='first').stack().astype(int)).mean()
    df_normalized=df.rank(method='min').stack().astype(int).map(rank_mean).unstack()
    return df_normalized
             
def col_row_norm(df):
    df = np.log2(df+1)
    df = df.subtract(df.median(),axis=1)
    df = df.subtract(df.mean(axis=1),axis=0)
    return df

def log_norm(df):
    df = np.log2(df+1)
    df = df.subtract(df.median(),axis=1)
    return df


def expr_with_info_plot(df,counts_row,outdir,xticklables,patient,metric='euclidean',method='average'):

    
    plt.figure(figsize = (3.5,5))
    gs = gridspec.GridSpec(1,1) 
    
    ## ==== heatmap of gene expression data
    loc=0  
    ax = plt.subplot(gs[loc,0])       
    g = sns.heatmap(df,cmap=plt.cm.PiYG_r,cbar=False,vmax=2,vmin=-2\
                    ,xticklabels=True,yticklabels=False,cbar_kws={"orientation": "vertical","use_gridspec":False}\
                    ,ax=ax)#"use_gridspec":False,"location":"top"})    
    # ax.xaxis.set_ticks_position('top')
    ax.set_xticklabels(xticklables,size = 16,rotation=45,ha='right')
    # ax.set_xlabel(None)
    ax.set_ylabel('\n\n{} genes'.format(df.shape[0]),fontsize=16)
    ax.yaxis.set_label_position('left')
    
    for i in np.arange(len(counts_row)+1):
        ax.axhline(y=sum(counts_row[:i]),xmin=-.02,xmax=1.02,color='k',lw=2.5,clip_on=False)    

    plt.title(patient,fontsize=16)     
#     plt.show()
    plt.savefig(outdir+'/{}.png'.format(patient),bbox_inches = 'tight',pad_inches=0.1,transparent=True,dpi=600)
    plt.close()





## ==== main section
    
outdir = 'f2_xenograft_expression_heatmap'
os.makedirs(outdir,exist_ok=True)

project_dir="/nv/vol190/zanglab/zw5j/wwork2018/AML_fran"
# project_dir="/Volumes/zanglab/zw5j/wwork2018/AML_fran"

# time serial data
expr_file='{}/manuscript_figs/f7_time_point_RNA/xenograft_new/f2_deg/f0_data_TPM/tpm_U4.csv'.format(project_dir)
expr_df = pd.read_csv(expr_file,index_col=0)

patient_datasets = {'U4':['XENOGRAPH7807','XENOGRAPH7811','XENOGRAPH7813','XENOGRAPH7791','XENOGRAPH7801','XENOGRAPH7794']}
patient_samples= ['U4']
xticklables = ['UT-rep1','UT-rep2','UT-rep3','T29-rep1','T29-rep2','T29-rep3']

# use PAML DEG
deg_clustering_dir = '{}/updated_202010/f1_kmeans_patients_by_deg/f1_kmeans_PAML_patients_by_PAML_DEG'.format(project_dir)


for patient in patient_samples:
    
    df = expr_df[patient_datasets[patient]]
    df = df[df.max(axis=1)>1].dropna()
    df = np.log2(df+1)
    deg_df = df
    
    # z-scored TPM
    df_zscored = pd.DataFrame(stats.zscore(df,axis=1))
    df_zscored.columns = df.columns
    df_zscored.index = df.index
    # print(df.shape)
    deg_df = df_zscored
    
    # log2 fold changes
    df_new = pd.DataFrame()
    for col in df.columns[1:]:
        df_new[col] = df[col]-df[df.columns[0]]
    # deg_df = df_new
    
    # PAML DEG
    deg_list1 = [i.strip() for i in open(deg_clustering_dir+'/txt/PAML_pthre5_rk3_ck2_genes_cluster1.txt').readlines()]
    deg_list2 = [i.strip() for i in open(deg_clustering_dir+'/txt/PAML_pthre5_rk3_ck2_genes_cluster2.txt').readlines()]
    deg_list3 = [i.strip() for i in open(deg_clustering_dir+'/txt/PAML_pthre5_rk3_ck2_genes_cluster3.txt').readlines()]
    # keep only those SG DEG
    deg_list1 = [i for i in deg_list1 if i in deg_df.index]
    deg_list2 = [i for i in deg_list2 if i in deg_df.index]
    deg_list3 = [i for i in deg_list3 if i in deg_df.index]
    
    # ==== re-order the genes followed by down-/div-/up-genes
    reindex = deg_list1+deg_list3+deg_list2
    # reindex = deg_list3
    deg_df = deg_df.reindex(reindex)
    
    # position for cluster box
    counts_row = [len(deg_list1),len(deg_list3),len(deg_list2)]
    # counts_row = [len(deg_list3)]
#     xticklables = timepoints[:df.shape[1]]
    expr_with_info_plot(deg_df,counts_row,outdir,xticklables,patient)
    
