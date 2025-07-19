import sys,argparse,re
import os,glob
import numpy as np
import pandas as pd
import operator

import matplotlib
from matplotlib import gridspec
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
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
from lifelines import KaplanMeierFitter
from lifelines import CoxPHFitter
from lifelines.utils import k_fold_cross_validation



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
        g1 = kmf_control.fit(t1, e1).plot(ax=ax,show_censors=True,label='SC low',\
                        censor_styles={'ms': 12, 'marker': '+'},ci_show=False,c='b',ls='-')
    
        kmf_exp = KaplanMeierFitter()
        #g2 = kmf_exp.fit(t2, e2, label=legends[1]).plot(ax=ax,show_censors=True,\
        g2 = kmf_exp.fit(t2, e2).plot(ax=ax,show_censors=True,label='SC high',\
                    censor_styles={'ms': 12, 'marker': '+'},ci_show=False,c='r',ls='-')
        handles, labels = ax.get_legend_handles_labels()
        # label_a_index = labels.index('SC score low')
        # label_b_index = labels.index('SC score high')
        # handles = np.append(handles[label_a_index],handles[label_b_index])
    #if pvalue<0.05:
        plt.axes().text(df['time'].max()*0.75,0.45,'p={:.2e}'.format(pvalue),fontsize=14,ha='center')
        plt.ylim([-0.02,1.05])
#     plt.xlim([0,max_val*1])
#         plt.title(title,fontsize=22)
        plt.xlabel('Days',fontsize=18)
        plt.ylabel('Survival probability',fontsize=18)
        plt.savefig(figname,bbox_inches='tight',pad_inches=.1,dpi=600,transparent=True)
        plt.show()
        plt.close()
    return results



# == main ==
outdir='f7b_cox_selected_div_TR_TCGA_AML_survival'
os.makedirs(outdir,exist_ok=True)


project_dir="/nv/vol190/zanglab/zw5j/wwork2018/AML_fran"
# project_dir="/Volumes/zanglab/zw5j/wwork2018/AML_fran"
# deg_clustering_dir = '{}/updated_202010/f1_kmeans_patients_by_deg/f1_kmeans_PAML_patients_by_PAML_DEG'.format(project_dir)
# patient_clustering_dir = '{}/updated_202010/f2_kmeans_patients_by_div_gene/f1_kmeans_PAML_patients_by_PAML_div_genes'.format(project_dir)
data_dir = '{}/f0_fran_data_new/data_processed/'.format(project_dir)
bart_results = '{}/updated_202010/f4_bart2_TR_expression/f2_PAML_DEG_DMG_removed_bart2/f1_DMG_removed_DEG_bart2_results/PAML_pthre5_rk3_ck2_genes_cluster3_div_bart_results.txt'.format(project_dir)
# beatAML_dir='{}/manuscript_figs/f11_clinical_TCGA_and_beatAML/beatAML/'.format(project_dir)
tcga_AML_dir='{}/manuscript_figs/f11_clinical_TCGA_and_beatAML/TCGA_LAML/'.format(project_dir)

# ==== TPM and clinical match ID
cohort = 'PAML' 
clinical_df = pd.read_csv(data_dir+os.sep+'{}_clinical.csv'.format(cohort),index_col=0)
tpm_df = pd.read_csv(data_dir+os.sep+'{}_TPM_diagnosis.txt'.format(cohort),index_col=0)

# ==== info of deg clustering
basename = '{}_pthre5_rk3_ck2'.format(cohort)
# deg_csv = deg_clustering_dir+'/csv/{}.csv'.format(basename)
# div_genes = [i.strip() for i in open(deg_clustering_dir+'/txt/{}_genes_cluster3.txt'.format(basename)).readlines()]
# top_tfs = pd.read_csv(bart_results,sep='\t',index_col=0).index[:20]
# top_tfs= set(['SPI1','STAG1','BCL6','CTCF','JMJD1C',
#          'GATA1','SMARCA4','RAD21','NFE2L1','TAL1','CTCF','SMC3',
#          'RBBP5','HDAC1','XRN2','GABPA'])
top_tfs= set(['GATA1','SMARCA4','RAD21','NFE2L1','TAL1','CTCF','SMC3'])
 

# prepare the matrix for cox regression
expression_data = np.transpose(tpm_df.loc[top_tfs])
cph_data = pd.concat([clinical_df['time_until_relapse_days'],expression_data],axis=1)
# cph_data = pd.concat([clinical_df,expression_data],axis=1)

# add penalizer term to the cox regression 
cph = CoxPHFitter()
# cph = CoxPHFitter(penalizer=0.1, l1_ratio=1.0) # sparse solutions,
cph.fit(cph_data, duration_col='time_until_relapse_days')
# print(cph.summary)
cph.summary.to_csv(outdir+os.sep+'cph_summary.csv')

# plot log HR
plt.figure(figsize=(4,4))
cph.plot()
plt.savefig(outdir+os.sep+'cph_logHR.pdf',bbox_inches='tight',pad_inches=.1,dpi=600,transparent=True)
plt.show()
plt.close()



# rank patients by TF20 score from beatAML expression


tcga_AML_clinical_df = pd.read_csv('{}/f1_clinical_of_target_TF/selected_case_id_with_clinical_info.csv'.format(tcga_AML_dir),index_col=0) 
tcga_AML_clinical_df = tcga_AML_clinical_df.fillna(0)
tcga_AML_clinical_df['time_max']=tcga_AML_clinical_df[['days_to_last_follow_up', 'days_to_death']].max(axis=1)

tcga_AML_expression_file='{}/f0_FPKM2TPM/TCGA_LAML_FPKM_Symbol.csv'.format(tcga_AML_dir)
tcga_AML_expression_df = pd.read_csv(tcga_AML_expression_file,index_col=0)
tcga_AML_expression_df = tcga_AML_expression_df[tcga_AML_clinical_df.index]

top_tfs = tcga_AML_expression_df.index.intersection(top_tfs)
tcga_AML_data=np.transpose(tcga_AML_expression_df.loc[top_tfs])

inner_product=np.inner(tcga_AML_data,cph.summary['coef'].loc[top_tfs])
inner_df = pd.DataFrame(inner_product)
inner_df.index = tcga_AML_data.index
c1_lower = inner_df.sort_values(by=[0],ascending=True).index[:int(.5*inner_df.shape[0])]
c2_higher = inner_df.sort_values(by=[0],ascending=True).index[-1*int(.5*inner_df.shape[0]):]


# plot the time-to-relapse
survival_df = pd.DataFrame(index = np.append(c1_lower,c2_higher))#;print(survival_df)
survival_df.loc[set(c1_lower).intersection(tcga_AML_clinical_df.index),'group']='c1'
survival_df.loc[set(c2_higher).intersection(tcga_AML_clinical_df.index),'group']='c2'
survival_df['status']= tcga_AML_clinical_df.loc[survival_df.index]['vital_status'] 
survival_df.loc[survival_df['status']=='Dead','status']=1
survival_df.loc[survival_df['status']=='Alive','status']=0

survival_df['time']= tcga_AML_clinical_df.loc[survival_df.index]['time_max']
survival_df = survival_df.dropna(axis=0,how='any')
avg_time_c1 = survival_df[survival_df['group']=='c1']['time'].mean()
avg_time_c2 = survival_df[survival_df['group']=='c2']['time'].mean()

figname=outdir+os.sep+'clinical_by_cox.png'
results = survival_for_two(survival_df,'c1','c2',figname)
survival_df.to_csv('{}/clinical_by_cox_info.csv'.format(outdir))
    
clinical_results = pd.DataFrame() 
clinical_results.loc[figname,'pvalue'] = results.p_value
clinical_results.loc[figname,'statistics'] = results.test_statistic
clinical_results.loc[figname,'c1_days'] = avg_time_c1
clinical_results.loc[figname,'c2_days'] = avg_time_c2

clinical_results.to_csv('{}/clinical_by_cox_LogRankSum.csv'.format(outdir))




# prepare the matrix for cox regression
cph_data = pd.concat([tcga_AML_clinical_df.loc[survival_df.index]['time_max'],tcga_AML_data.loc[survival_df.index]],axis=1)
# cph_data = pd.concat([clinical_df,expression_data],axis=1)

# add penalizer term to the cox regression 
cph = CoxPHFitter()
# cph = CoxPHFitter(penalizer=0.1, l1_ratio=1.0) # sparse solutions,
cph.fit(cph_data, duration_col='time_max')
# print(cph.summary)
cph.summary.to_csv(outdir+os.sep+'cph_summary_tcga_AML.csv')

# plot log HR
plt.figure(figsize=(4,4))
cph.plot()
plt.savefig(outdir+os.sep+'cph_logHR_tcga_AML.pdf',bbox_inches='tight',pad_inches=.1,dpi=600,transparent=True)
plt.show()
plt.close()



