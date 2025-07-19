import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
#from GenomeData import *
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=10
import seaborn as sns
sns.set(font_scale=1)
sns.set_style("whitegrid", {'axes.grid' : False})
# sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
from matplotlib import gridspec
#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
from  scipy.cluster.hierarchy import fcluster    



def return_patient_info(project_dir,cohort):
    
    '''This is used to return the deseq2 results and TPM for patient samples '''
    file_dir = '{}/f0_fran_data_new/data_processed'.format(project_dir)
    
    tpm_diagnosis = pd.read_csv(file_dir+os.sep+'{}_TPM_diagnosis.txt'.format(cohort),sep=',',index_col=0)
    tpm_relapse = pd.read_csv(file_dir+os.sep+'{}_TPM_relapse.txt'.format(cohort),sep=',',index_col=0)
    tpm_normal = pd.read_csv(file_dir+os.sep+'{}_TPM_normal.txt'.format(cohort),sep=',',index_col=0)
    log2FC = pd.read_csv(file_dir+os.sep+'{}_deseq2_log2FC.txt'.format(cohort),sep=',',index_col=0)
    
    # return log2FC, np.log2(tpm_diagnosis+1), np.log2(tpm_relapse+1)
    return log2FC, tpm_diagnosis, tpm_relapse


## ==== main section
    
# outdir = 'f1_PAML_heatmap_with_mutation'
# os.makedirs(outdir,exist_ok=True)

# project_dir="/nv/vol190/zanglab/zw5j/wwork2018/AML_fran"
project_dir="/Volumes/zanglab/zw5j/wwork2018/AML_fran"
clinical_dir = '{}/f0_fran_data_new/data_processed/'.format(project_dir)

cohort = 'PAML'
gene_list = ['DNMT3B',
 'ZBTB46',
 'NYNRIN',
 'ARHGAP22',
 'LAPTM4B',
 'MMRN1',
 'DPYSL3',
 'KIAA0125',
 'CDK6',
 'CPXM1',
 'SOCS2',
 'C19orf77',
 'EMP1',
 'NGFRAP1',
 'CD34',
 'AKR1C3',
 'GPR56']

# ==== gene expression TPM
log2FC, tpm_diagnosis, tpm_relapse =  return_patient_info(project_dir,cohort)
info = pd.read_csv(clinical_dir+os.sep+'{}_clinical.csv'.format(cohort),index_col=0)
df_sub = log2FC.loc[gene_list]
df_sub.columns = info.loc[df_sub.columns]['subject_ID']
percentage_each_patient = 100*(df_sub>0).sum(axis=1).sort_values(ascending=False)/df_sub.shape[1]
df_percentage = pd.DataFrame(percentage_each_patient,columns=['% up patients']).round(2)
pd.concat([df_sub.round(2),df_percentage],axis=1).to_csv('log2FC_17_LSC_signature_genes.csv')

   
    
    
