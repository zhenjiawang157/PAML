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



def survival_for_two(df,treat,ctrl,figname,thre):
    
    # select the time and status info for treat and control group
    ix = df['group'] == treat
    
    t1 = df.loc[ix]['time']
    e1 = df.loc[ix]['status'] 
    t2 = df.loc[~ix]['time']
    e2 = df.loc[~ix]['status']
    
    results = logrank_test(t1,t2,event_observed_A = e1,event_observed_B = e2)
    pvalue = results.p_value#;print('pvalue:\t{}'.format(pvalue))

    # survival curves
    plt.figure(figsize=(4.,3))
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
    plt.xticks([0,250,500,750,1000,1250,1500])
#     plt.xlim([0,max_val*1])
    # plt.title('TR pthre {}'.format(np.power(.1,p_thre).round(3)),fontsize=16)
    plt.xlabel('Days',fontsize=16)
    plt.ylabel('Survival probability',fontsize=16)
    plt.legend(bbox_to_anchor=[1.0,1.0],loc='upper right',frameon=False,\
               borderaxespad=0.2,labelspacing=.2,handletextpad=0.2,handlelength=1)
    plt.savefig(figname,bbox_inches='tight',pad_inches=.1,dpi=300,transparent=True)
    plt.show()
    plt.close()
    
    return results


def pre_data(project_dir):
    # ==== beat AML data
    beatAML_dir='{}/manuscript_figs/f11_clinical_TCGA_and_beatAML/beatAML/'.format(project_dir)
    beatAML_expression_file='{}/f0_expression_TPM/beatAML_RPKM_Symbol.csv'.format(beatAML_dir)
    beatAML_expression_df = pd.read_csv(beatAML_expression_file,index_col=0)
    beatAML_clinical_df_ori = pd.read_excel('{}/beatAML_suppl_ClinicalInformation_7+3Tx.xlsx'.format(beatAML_dir),index_col=0,sheet_name='7+3') 
    beatAML_clinical_df = beatAML_clinical_df_ori[['vitalStatus','overallSurvival']].dropna()
    
    # ==== check the patients with both clinical data and expression data
    shared_patients = beatAML_clinical_df.index.intersection(beatAML_expression_df.columns)
    beatAML_clinical_df_ori.loc[shared_patients].to_csv(outdir+os.sep+'beatAML_suppl_ClinicalInformation_7+3Tx_new.csv')
    beatAML_expression_df = beatAML_expression_df[shared_patients]
    beatAML_clinical_df = beatAML_clinical_df.loc[shared_patients]
    beatAML_clinical_df.loc[beatAML_clinical_df['vitalStatus']=='Dead','vitalStatus']=1
    beatAML_clinical_df.loc[beatAML_clinical_df['vitalStatus']=='Alive','vitalStatus']=0
    
    return beatAML_expression_df, beatAML_clinical_df



# == main ==
outdir='f4_cox_on_pvalueThre_div_TR_beatAML_survival_CV'
os.makedirs(outdir,exist_ok=True)

project_dir="/nv/vol190/zanglab/zw5j/wwork2018/AML_fran"
project_dir="/Volumes/zanglab/zw5j/wwork2018/AML_fran"
beatAML_expression_df, beatAML_clinical_df = pre_data(project_dir)

# ==== get the BART results for divergent genes
bart_results = '{}/updated_202010/f4_bart2_TR_expression/f2_PAML_DEG_DMG_removed_bart2/f1_DMG_removed_DEG_bart2_results/PAML_pthre5_rk3_ck2_genes_cluster3_div_bart_results.txt'.format(project_dir)
tfs_df = pd.read_csv(bart_results,sep='\t',index_col=0)
p_thre = 3
top_tfs = tfs_df[tfs_df['irwin_hall_pvalue']<np.power(.1,p_thre)].index;print(len(top_tfs))
# print(len(top_tfs));continue

# ==== repare the matrix for cox regression
beatAML_data=np.transpose(beatAML_expression_df.loc[top_tfs].dropna()) # gene on the column
cph_data = pd.concat([beatAML_clinical_df,beatAML_data],axis=1)
penalizers = [.001, .002, .005, .01, .02, .05, .1, .2, .5]
l1_ratios = [0, .2, .4, .6, .8, 1]
summary_score=pd.DataFrame()

for penalizer in penalizers[:1]:
    for l1_ratio in l1_ratios[2:3]:
        basename = 'TR_Pthre{}_Cox_penalizer{}_l1_ratios{}'.format(p_thre,penalizer,l1_ratio)
        print(basename)
        ## ==== cross valication to test the parameters
        k = 5
        cph = CoxPHFitter(penalizer=penalizer, l1_ratio=l1_ratio) # sparse solutions,
        for ii in np.arange(10):
            scores = k_fold_cross_validation(cph, cph_data,'overallSurvival',event_col='vitalStatus',k=k,scoring_method="concordance_index")
            summary_score.loc[basename,'cv_time{}_mean'.format(ii)]=np.mean(scores).round(4)
        cv_cols = [i for i in summary_score.columns if re.search('cv_time',i)]
        summary_score.loc[basename,'cv_all_mean'] = summary_score.loc[basename,cv_cols].mean()
        summary_score.loc[basename,'cv_all_max'] = summary_score.loc[basename,cv_cols].max()

        # cph = CoxPHFitter()
        cph = CoxPHFitter(penalizer=penalizer, l1_ratio=l1_ratio) # sparse solutions,
        cph.fit(cph_data, duration_col='overallSurvival', event_col='vitalStatus')
        cph.summary.to_csv(outdir+os.sep+'{}_cph_summary.csv'.format(basename))
        # print(cph.summary)
        # ==== plot log HR
        selected_genes = cph.summary['coef'][np.abs(cph.summary['coef'])>0].sort_values(ascending=False).index
        plt.figure(figsize=(3,7))
        # cph.plot()
        cph.plot(selected_genes)
        plt.ylim([-1,len(selected_genes)])
        # plt.title('TR pthre {}'.format(np.power(.1,p_thre).round(3)),fontsize=16)
        plt.savefig(outdir+os.sep+'{}_logHR.pdf'.format(basename),bbox_inches='tight',pad_inches=.1,dpi=600,transparent=True)
        plt.show()
        plt.close()
        # print(len(selected_genes))


        ## ==== survival anlaysis
        inner_product=np.inner(beatAML_data,cph.summary['coef'])
        inner_df = pd.DataFrame(inner_product)
        inner_df.index = beatAML_data.index
        c1_lower = inner_df.sort_values(by=[0],ascending=True).index[:int(.5*inner_df.shape[0])]
        c2_higher = inner_df.sort_values(by=[0],ascending=True).index[-1*int(.5*inner_df.shape[0]):]
        # plot the clinical figs
        survival_df = pd.DataFrame(index = np.append(c1_lower,c2_higher))#;print(survival_df)
        survival_df.loc[set(c1_lower).intersection(beatAML_clinical_df.index),'group']='c1'
        survival_df.loc[set(c2_higher).intersection(beatAML_clinical_df.index),'group']='c2'
        survival_df['status']= beatAML_clinical_df.loc[survival_df.index]['vitalStatus']         
        survival_df['time']= beatAML_clinical_df.loc[survival_df.index]['overallSurvival']
        survival_df = survival_df.dropna(axis=0,how='any')
        avg_time_c1 = survival_df[survival_df['group']=='c1']['time'].mean()
        avg_time_c2 = survival_df[survival_df['group']=='c2']['time'].mean()
        figname=outdir+os.sep+'{}_clinical.pdf'.format(basename)
        results = survival_for_two(survival_df,'c1','c2',figname,p_thre)
        survival_df.to_csv('{}/{}_clinical_info.csv'.format(outdir,basename))
        # data for clinical analysis
        summary_score.loc[basename,'logrank_pvalue'] = results.p_value
        summary_score.loc[basename,'logrank_statistics'] = results.test_statistic
        summary_score.loc[basename,'logrank_c1_low_days'] = avg_time_c1
        summary_score.loc[basename,'logrank_c2_high_days'] = avg_time_c2
        ## ==== life table
        life_table = pd.DataFrame()
        risk_days = [0,250,500,750,1000,1250,1500]
        for risk_day in risk_days:
            life_table.loc['Bottom 50%',risk_day] = survival_df[(survival_df['group']=='c1') & (survival_df['time']>=risk_day)].shape[0]
            life_table.loc['Top 50%',risk_day] = survival_df[(survival_df['group']=='c2') & (survival_df['time']>=risk_day)].shape[0]
        life_table.to_csv('{}/{}_life_table.csv'.format(outdir,basename))
        
## == save out the results
summary_score.to_csv('{}/summary_score.csv'.format(outdir))        



## ==== adaptive lasso
# x = beatAML_data
# gprime = lambda w: 1.0/(2.*np.sqrt(np.abs(w))+np.finfo(float).eps)
# n_samples,n_features = x.shape
# weights = np.ones(n_features)
# selected_features = n_features
# for ii in np.arange(10):
#     X_w = x / weights
#     cph_data = pd.concat([beatAML_clinical_df,X_w],axis=1)
#     cph = CoxPHFitter(penalizer=.01, l1_ratio=1) # sparse solutions,
#     cph.fit(cph_data, duration_col='overallSurvival', event_col='vitalStatus')
#     # print(cph.summary)
#     coef_ = cph.summary['coef']/weights
#     weights = gprime(coef_)
#     selected_features = len([i for i in coef_ if i !=0])
#     # cph.summary.to_csv(outdir+os.sep+'Pthre{}_cph_summary_beatAML.csv'.format(p_thre))
#     print('{}, {}'.format(ii,selected_features))
#     # ==== plot log HR
#     # selected_genes = cph.summary['coef'][np.abs(cph.summary['coef'])>0].sort_values(ascending=False).index
#     plt.figure(figsize=(3,7))
#     cph.plot()
#     # cph.plot(selected_genes)
#     plt.title('TR pthre {}'.format(np.power(.1,p_thre).round(3)),fontsize=16)
#     # plt.savefig(outdir+os.sep+'Pthre{}_cph_logHR_beatAML.pdf'.format(p_thre),bbox_inches='tight',pad_inches=.1,dpi=600,transparent=True)
#     plt.show()
#     plt.close()
    
