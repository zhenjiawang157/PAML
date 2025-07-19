import sys,argparse,re
import os,glob
import numpy as np
import pandas as pd
import operator

import matplotlib
from matplotlib import gridspec
# matplotlib.use('Agg')
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

    # survival curves
    plt.figure(figsize=(3.,3.))
    ax = plt.subplot(111)
    kmf_exp = KaplanMeierFitter()
    #g2 = kmf_exp.fit(t2, e2, label=legends[1]).plot(ax=ax,show_censors=True,\
    g2 = kmf_exp.fit(t2, e2).plot(ax=ax,show_censors=True,label='Top 50%',\
                censor_styles={'ms': 12, 'marker': '+'},ci_show=False,c='r',ls='-')
    
    kmf_control = KaplanMeierFitter()
    #g1 = kmf_control.fit(t1, e1, label=legends[0]).plot(ax=ax,show_censors=True,\
    g1 = kmf_control.fit(t1, e1).plot(ax=ax,show_censors=True,label='Bottom 50%',\
                    censor_styles={'ms': 12, 'marker': '+'},ci_show=False,c='b',ls='-')

    plt.axes().text(df['time'].max()*0.75,0.45,'p={:.2e}'.format(pvalue),fontsize=12,ha='center')
    plt.ylim([-0.02,1.05])
#     plt.xlim([0,max_val*1])
#         plt.title(title,fontsize=22)
    plt.xlabel('Days',fontsize=16)
    plt.ylabel('Survival probability',fontsize=16)
    plt.legend(bbox_to_anchor=[1.03,1.03],loc='upper right',frameon=False,\
               borderaxespad=0.2,labelspacing=.2,handletextpad=0.2,handlelength=1)
    plt.savefig(figname,bbox_inches='tight',pad_inches=.1,dpi=300,transparent=True)
    plt.show()
    plt.close()
    
    return results



# == main ==
outdir='f1_cox_on_top_div_TR_beatAML_survival'
os.makedirs(outdir,exist_ok=True)

project_dir="/nv/vol190/zanglab/zw5j/wwork2018/AML_fran"
project_dir="/Volumes/zanglab/zw5j/wwork2018/AML_fran"
# deg_clustering_dir = '{}/updated_202010/f1_kmeans_patients_by_deg/f1_kmeans_PAML_patients_by_PAML_DEG'.format(project_dir)
# div_genes = [i.strip() for i in open(deg_clustering_dir+'/txt/PAML_pthre5_rk3_ck2_genes_cluster3.txt').readlines()]
# data_dir = '{}/f0_fran_data_new/data_processed/'.format(project_dir)

# ==== beat AML data
beatAML_dir='{}/manuscript_figs/f11_clinical_TCGA_and_beatAML/beatAML/'.format(project_dir)
beatAML_expression_file='{}/f0_expression_TPM/beatAML_RPKM_Symbol.csv'.format(beatAML_dir)
beatAML_expression_df = pd.read_csv(beatAML_expression_file,index_col=0)
beatAML_clinical_df = pd.read_excel('{}/beatAML_suppl_ClinicalInformation_7+3Tx.xlsx'.format(beatAML_dir),index_col=0,sheet_name='7+3') 
beatAML_clinical_df = beatAML_clinical_df[['vitalStatus','overallSurvival']].dropna()
# ==== check the patients with both clinical data and expression data
shared_patients = beatAML_clinical_df.index.intersection(beatAML_expression_df.columns)
beatAML_expression_df = beatAML_expression_df[shared_patients]
beatAML_clinical_df = beatAML_clinical_df.loc[shared_patients]
beatAML_clinical_df.loc[beatAML_clinical_df['vitalStatus']=='Dead','vitalStatus']=1
beatAML_clinical_df.loc[beatAML_clinical_df['vitalStatus']=='Alive','vitalStatus']=0

# ==== get the BART results for divergent genes
bart_results = '{}/updated_202010/f4_bart2_TR_expression/f2_PAML_DEG_DMG_removed_bart2/f1_DMG_removed_DEG_bart2_results/PAML_pthre5_rk3_ck2_genes_cluster3_div_bart_results.txt'.format(project_dir)
tfs_df = pd.read_csv(bart_results,sep='\t',index_col=0)
# top_tfs= set(['GATA1','SMARCA4','RAD21','NFE2L1','TAL1','CTCF','SMC3'])
# top_tfs= set(['JMJD1C','TEAD4','ERG','ZBTB7A','ETS1'])
top_thres=[5,10,15,20]

for top_thre in top_thres[:]:
    # ==== info of deg clustering
    top_tfs = tfs_df.index[:top_thre]
    beatAML_data=np.transpose(beatAML_expression_df.loc[top_tfs].dropna()) # gene on the column

    # p==== repare the matrix for cox regression
    cph_data = pd.concat([beatAML_clinical_df,beatAML_data],axis=1)
    cph = CoxPHFitter()
    # cph = CoxPHFitter(penalizer=0.2, l1_ratio=1) # sparse solutions,
    cph.fit(cph_data, duration_col='overallSurvival', event_col='vitalStatus')
    # print(cph.summary)
    cph.summary.to_csv(outdir+os.sep+'TRthre{}_cph_summary_beatAML.csv'.format(top_thre))
    # selected_genes = cph.summary['coef'][np.abs(cph.summary['coef'])>0.01].index
    # plot log HR
    plt.figure(figsize=(4,4))
    cph.plot()
    plt.savefig(outdir+os.sep+'TRthre{}_cph_logHR_beatAML.pdf'.format(top_thre),bbox_inches='tight',pad_inches=.1,dpi=600,transparent=True)
    plt.show()
    plt.close()
    print(cph.summary['coef'])

    # k=5
    # cv_score=pd.DataFrame(columns=np.arange(k))
    # cph_cv = CoxPHFitter()
    # for ii in np.arange(10):
    #     scores = k_fold_cross_validation(cph_cv, cph_data,'overallSurvival',event_col='vitalStatus',k=k,scoring_method="concordance_index")
    #     # scores = k_fold_cross_validation(cph, cph_data[np.append(selected_genes,cph_data.columns[:2])],\
    #     #                               'overallSurvival', event_col='vitalStatus', k=k,scoring_method="concordance_index")
    #     cv_score.loc[ii]=scores
    # cv_score.to_csv(outdir+os.sep+'CV_score.csv')
    # print(cph.summary['coef'])
    


    # ==== survival anlaysis
    inner_product=np.inner(beatAML_data,cph.summary['coef'])
    inner_df = pd.DataFrame(inner_product)
    inner_df.index = beatAML_data.index
    c1_lower = inner_df.sort_values(by=[0],ascending=True).index[:int(.5*inner_df.shape[0])]
    c2_higher = inner_df.sort_values(by=[0],ascending=True).index[-1*int(.5*inner_df.shape[0]):]
    
    # plot the time-to-relapse
    survival_df = pd.DataFrame(index = np.append(c1_lower,c2_higher))#;print(survival_df)
    survival_df.loc[set(c1_lower).intersection(beatAML_clinical_df.index),'group']='c1'
    survival_df.loc[set(c2_higher).intersection(beatAML_clinical_df.index),'group']='c2'
    survival_df['status']= beatAML_clinical_df.loc[survival_df.index]['vitalStatus'] 
    # survival_df.loc[survival_df['status']=='Dead','status']=1
    # survival_df.loc[survival_df['status']=='Alive','status']=0
    
    survival_df['time']= beatAML_clinical_df.loc[survival_df.index]['overallSurvival']
    survival_df = survival_df.dropna(axis=0,how='any')
    avg_time_c1 = survival_df[survival_df['group']=='c1']['time'].mean()
    avg_time_c2 = survival_df[survival_df['group']=='c2']['time'].mean()
    
    figname=outdir+os.sep+'TRthre{}_clinical_by_cox.png'.format(top_thre)
    results = survival_for_two(survival_df,'c1','c2',figname)
    survival_df.to_csv('{}/TRthre{}_clinical_by_cox_info.csv'.format(outdir,top_thre))
    
# clinical_results = pd.DataFrame() 
# clinical_results.loc[figname,'pvalue'] = results.p_value
# clinical_results.loc[figname,'statistics'] = results.test_statistic
# clinical_results.loc[figname,'c1_days'] = avg_time_c1
# clinical_results.loc[figname,'c2_days'] = avg_time_c2
# clinical_results.to_csv('{}/clinical_by_cox_LogRankSum.csv'.format(outdir))


