import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
from scipy import stats
#from GenomeData import *
import matplotlib
# matplotlib.use('Agg')
from matplotlib import gridspec
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=18
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")
#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')
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



def plot_fillbetween(patient,ax,positions,deg_df,current_genes,color,ylabel):
    
    # ==== for each if down/div/up genes, plot their expression pattern
    max_values,top_values,middle_values,bottom_values,min_values = [],[],[],[],[]
    for col in deg_df.columns:
        expr_values = deg_df.loc[current_genes][col].values
        max_values.append(np.percentile(expr_values,90))
        top_values.append(np.percentile(expr_values,75))
        middle_values.append(np.percentile(expr_values,50))
        bottom_values.append(np.percentile(expr_values,25))
        min_values.append(np.percentile(expr_values,10))
        
    ax.scatter(positions,middle_values,c=color)
    ax.plot(middle_values,c=color)
    # ax.fill_between(positions,max_values,min_values,color=color,interpolate=True,lw=0,alpha=.2)
    ax.fill_between(positions,top_values,bottom_values,color=color,interpolate=True,lw=0,alpha=.5)
    ax.set_ylabel(ylabel,linespacing = 2)
    if patient =='AML130':
        ax.set_ylim([0.68,4.9])
    if patient =='AML124':
        ax.set_ylim([0.8,4.4])
    ax.set_yticks([1,2,3,4])
    # ax.set_ylim([0.,5.85])




## ==== main section
    
outdir = 'f1_time_serial_expression_trend_QN_log2TPM'
os.makedirs(outdir,exist_ok=True)

project_dir="/nv/vol190/zanglab/zw5j/wwork2018/AML_fran"
project_dir="/Volumes/zanglab/zw5j/wwork2018/AML_fran"

# time serial data
expr_file='{}/manuscript_figs/f7_time_point_RNA/bulk_serial_new//f0_process/f3_expr/tpm.csv'.format(project_dir)
expr_df = pd.read_csv(expr_file,index_col=0)

patient_datasets = {'AML130':['1245','1477','1575','1819','1926'],'AML124':['882','1482','1594'],'AML126':['969','1435','1598','1730']}
timepoints = ['Diagnosis','Time point1','Time point2','Time point3','Time point4']
patient_samples=['AML130','AML124','AML126']

# use PAML DEG
deg_clustering_dir = '{}/updated_202010/f1_kmeans_patients_by_deg/f1_kmeans_PAML_patients_by_PAML_DEG'.format(project_dir)


for patient in patient_samples[:]:
    # expression values
    df = expr_df[patient_datasets[patient]]
    df = df[df.max(axis=1)>1].dropna()
    df = np.log2(df+1)
    # deg_df = df
    deg_df = quantile_normalization(df)
    print(patient,deg_df.shape)  
      
    # PAML DEG
    deg_list1 = [i.strip() for i in open(deg_clustering_dir+'/txt/PAML_pthre5_rk3_ck2_genes_cluster1.txt').readlines()]
    deg_list2 = [i.strip() for i in open(deg_clustering_dir+'/txt/PAML_pthre5_rk3_ck2_genes_cluster2.txt').readlines()]
    deg_list3 = [i.strip() for i in open(deg_clustering_dir+'/txt/PAML_pthre5_rk3_ck2_genes_cluster3.txt').readlines()]
    # keep only those SG DEG
    deg_list1 = [i for i in deg_list1 if i in deg_df.index]
    deg_list2 = [i for i in deg_list2 if i in deg_df.index]
    deg_list3 = [i for i in deg_list3 if i in deg_df.index]
    
    # ==== re-order the genes followed by down-/div-/up-genes
    # reindex = deg_list1+deg_list3+deg_list2
    # deg_df = deg_df.reindex(reindex)
        
    # ==== plot the fig    
    plt.figure()
    f, (ax1, ax2, ax3) = plt.subplots(3,sharex=True,figsize = (2.5,4))
    f.subplots_adjust(hspace = 0) 
    
    positions = np.arange(df.shape[1])
    xticklables = timepoints[:df.shape[1]]
    plot_fillbetween(patient,ax1,positions,deg_df,deg_list1,'tab:green','Group A')
    plot_fillbetween(patient,ax2,positions,deg_df,deg_list3,'tab:purple','Normalized expression \n Group B')
    plot_fillbetween(patient,ax3,positions,deg_df,deg_list2,'tab:red','Group C')
    
    ax3.set_xticks(positions)
    ax3.set_xticklabels(xticklables,rotation=45,ha='center',fontsize=15)
    
    ax1.set_title(patient)  
    plt.savefig(outdir+'/{}.pdf'.format(patient),bbox_inches = 'tight',pad_inches=0.1,transparent=True,dpi=600)
    plt.show()
    plt.close()

