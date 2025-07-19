'''
This file is used to plot the heatmap to compare the gene expression change from N to D to R,
by separeting the PAML/SG patients, which have RNA-seq data from different batches 
'''

import sys,argparse,re
import os,glob
import numpy as np
import pandas as pd
import operator
from collections import Counter
#from GenomeData import *
#import site_align_heatmap
from sklearn.cluster import KMeans
from sklearn.metrics.cluster import adjusted_rand_score

import matplotlib
from matplotlib import gridspec
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=14
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")
#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')
from lifelines.statistics import logrank_test
#from lifelines.plotting import add_at_risk_counts
from lifelines import KaplanMeierFitter


def survival_for_two(df,treat,ctrl,figname):
    
    # select the time and status info for treat and control group
    ix = df['group'] == treat
    
    t1 = df.loc[ix]['time']
    e1 = df.loc[ix]['status'] 
    t2 = df.loc[~ix]['time']
    e2 = df.loc[~ix]['status']
    
    results = logrank_test(t1,t2,event_observed_A = e1,event_observed_B = e2)
    pvalue = results.p_value#;print('pvalue:\t{}'.format(pvalue))

#    if pvalue<1e-7:
    if 1:
    # survival curves
        plt.figure(figsize=(3.,3.))
        ax = plt.subplot(111)
        kmf_control = KaplanMeierFitter()
        #g1 = kmf_control.fit(t1, e1, label=legends[0]).plot(ax=ax,show_censors=True,\
        g1 = kmf_control.fit(t1, e1).plot(ax=ax,show_censors=True,label='Patient C1',\
                        censor_styles={'ms': 12, 'marker': '+'},ci_show=False,c='tab:orange',ls='-')
    
        kmf_exp = KaplanMeierFitter()
        #g2 = kmf_exp.fit(t2, e2, label=legends[1]).plot(ax=ax,show_censors=True,\
        g2 = kmf_exp.fit(t2, e2).plot(ax=ax,show_censors=True,label='Patient C2',\
                    censor_styles={'ms': 12, 'marker': '+'},ci_show=False,c='tab:blue',ls='-')
        handles, labels = ax.get_legend_handles_labels()
        label_a_index = labels.index('Patient C1')
        label_b_index = labels.index('Patient C2')
        handles = np.append(handles[label_a_index],handles[label_b_index])
    #;print(labels)
#         lg = ax.legend(handles, legends,loc=0,fontsize=12,borderaxespad=0,handletextpad=.2,labelspacing=.1,handlelength=1.1,frameon=False)
    #if pvalue<0.05:
        plt.axes().text(df['time'].max()*0.75,0.45,'p={:.2f}'.format(pvalue),fontsize=14,ha='center')
        plt.ylim([-0.02,1.05])
#     plt.xlim([0,max_val*1])
#         plt.title(title,fontsize=22)
        plt.xlabel('Time to relapse (days)',fontsize=18)
        plt.ylabel('Relapse probability',fontsize=18)
        plt.savefig(figname,bbox_inches='tight',pad_inches=.1,dpi=600,transparent=True)
        plt.show()
#        plt.close()
    return results


def return_patient_info(cohort,logfc_thre=1,qvalue_thre=0.01):
    
    '''This is used to return the deseq2 results and TPM for patient samples '''
    file_dir = '/nv/vol190/zanglab/zw5j/wwork2018/AML_fran/f0_fran_data_new/data_processed'
#     file_dir = '/Volumes/zanglab/zw5j/wwork2018/AML_fran/f0_fran_data_new/data_processed'
    clinical_df = pd.read_csv(file_dir+os.sep+'{}_clinical.csv'.format(cohort),sep=',',index_col=0)
    
    return clinical_df



        
    
indirs = ['f1_kmeans_PAML_patients_by_PAML_div_genes','f2_kmeans_SG_patients_by_PAML_div_genes','f3_kmeans_SG_patients_by_PAML_div_SG_DEG']
indir_cohorts=['PAML','SG','SG']
outdir = 'f5_patient_clinical_merge_cohorts'
os.makedirs(outdir,exist_ok=True)
    
fig_basename = 'cohorts_merged'
clinical_results = pd.DataFrame()   
    
# PAML data
cohort='PAML'
PAML_clinical_df = return_patient_info(cohort)
PAML_patient_c1 = [i.strip() for i in open(indirs[0]+'/txt/{}_pthre5_rk3_ck2_patients_cluster1.txt'.format(cohort)).readlines()]
PAML_patient_c2 = [i.strip() for i in open(indirs[0]+'/txt/{}_pthre5_rk3_ck2_patients_cluster2.txt'.format(cohort)).readlines()]

# SG data
cohort='SG'
SG_clinical_df = return_patient_info(cohort)
SG_patient_c1 = [i.strip() for i in open(indirs[2]+'/txt/{}_pthre5_rk3_ck2_patients_cluster1.txt'.format(cohort)).readlines()]
SG_patient_c2 = [i.strip() for i in open(indirs[2]+'/txt/{}_pthre5_rk3_ck2_patients_cluster2.txt'.format(cohort)).readlines()]

# merged data
clinical_df = pd.concat([PAML_clinical_df,SG_clinical_df])
patient_c1 = PAML_patient_c2+SG_patient_c2
patient_c2 = PAML_patient_c1+SG_patient_c1


# clinical analysis
survival_df = pd.DataFrame(index = np.append(patient_c1,patient_c2))#;print(survival_df)
survival_df.loc[set(patient_c1).intersection(clinical_df.index),'group']='c1'
survival_df.loc[set(patient_c2).intersection(clinical_df.index),'group']='c2'
survival_df['status']= clinical_df.loc[survival_df.index]['age_at_diagnosis'] 
survival_df['time']= clinical_df.loc[survival_df.index]['time_until_relapse_days']
survival_df = survival_df.dropna(axis=0,how='any')
figname=outdir+os.sep+'{}.png'.format(fig_basename)
results = survival_for_two(survival_df,'c1','c2',figname)
survival_df.to_csv('{}/{}_survival_info.csv'.format(outdir,fig_basename))

clinical_results.loc[figname,'pvalue'] = results.p_value
clinical_results.loc[figname,'statistics'] = results.test_statistic
    
clinical_results.to_csv('{}/LogRankSum_survival_results.csv'.format(outdir))






        


