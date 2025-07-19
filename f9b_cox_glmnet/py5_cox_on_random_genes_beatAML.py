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
sns.set(font_scale=1.1)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")
#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')
from lifelines.statistics import logrank_test
from lifelines import KaplanMeierFitter
from lifelines import CoxPHFitter
from lifelines.utils import k_fold_cross_validation



def survival_for_two(df,treat,ctrl,figname,risk_days,flag):
    
    # select the time and status info for treat and control group
    ix = df['group'] == treat
    
    t1 = df.loc[ix]['time']
    e1 = df.loc[ix]['status'] 
    t2 = df.loc[~ix]['time']
    e2 = df.loc[~ix]['status']
    
    results = logrank_test(t1,t2,event_observed_A = e1,event_observed_B = e2)
    pvalue = results.p_value#;print('pvalue:\t{}'.format(pvalue))

    # survival curves
    plt.figure(figsize=(3.6,3.))
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
    plt.title('{}\n{}'.format(basename.split('Cox_')[1],flag),fontsize=15)
    plt.xlabel('Days',fontsize=15)
    plt.ylabel('Survival probability',fontsize=15)
    plt.xticks(risk_days)
    plt.legend(bbox_to_anchor=[1.0,1.0],loc='upper right',frameon=False,\
               borderaxespad=0.2,labelspacing=.2,handletextpad=0.2,handlelength=1)
    plt.savefig(figname,bbox_inches='tight',pad_inches=.1,dpi=300,transparent=True)
    plt.show()
    plt.close()
    
    return results


def calinical_analysis(summary_score,cph,logrank_expr_data,logrank_clinical_data,logrank_clinical_cols,flag,basename,risk_days):
    ## ==== survival anlaysis
    inner_product=np.inner(logrank_expr_data,cph.summary['coef'].loc[logrank_expr_data.columns])
    inner_df = pd.DataFrame(inner_product)
    inner_df.index = logrank_expr_data.index
    c1_lower = inner_df.sort_values(by=[0],ascending=True).index[:int(.5*inner_df.shape[0])]
    c2_higher = inner_df.sort_values(by=[0],ascending=True).index[-1*int(.5*inner_df.shape[0]):]
    # plot the clinical figs
    survival_df = pd.DataFrame(index = np.append(c1_lower,c2_higher))#;print(survival_df)
    survival_df.loc[set(c1_lower).intersection(logrank_clinical_data.index),'group']='c1'
    survival_df.loc[set(c2_higher).intersection(logrank_clinical_data.index),'group']='c2'
    survival_df['status']= logrank_clinical_data.loc[survival_df.index][logrank_clinical_cols[0]]         
    survival_df['time']= logrank_clinical_data.loc[survival_df.index][logrank_clinical_cols[1]]
    survival_df = survival_df.dropna(axis=0,how='any')
    avg_time_c1 = survival_df[survival_df['group']=='c1']['time'].mean()
    avg_time_c2 = survival_df[survival_df['group']=='c2']['time'].mean()
    figname=outdir+os.sep+'{}_clinical_{}.pdf'.format(basename,flag)
    results = survival_for_two(survival_df,'c1','c2',figname,risk_days,flag)
    survival_df.to_csv('{}/{}_clinical_{}_data.csv'.format(outdir,basename,flag))
    # data for clinical analysis
    summary_score.loc[basename,'{}_logrank_pvalue'.format(flag)] = results.p_value
    summary_score.loc[basename,'{}_logrank_statistics'.format(flag)] = results.test_statistic
    summary_score.loc[basename,'{}_logrank_c1_low_days'.format(flag)] = avg_time_c1
    summary_score.loc[basename,'{}_logrank_c2_high_days'.format(flag)] = avg_time_c2
    ## ==== life table 
    life_table = pd.DataFrame()
    for risk_day in risk_days:
        life_table.loc['Bottom 50%',risk_day] = survival_df[(survival_df['group']=='c1') & (survival_df['time']>=risk_day)].shape[0]
        life_table.loc['Top 50%',risk_day] = survival_df[(survival_df['group']=='c2') & (survival_df['time']>=risk_day)].shape[0]
    life_table.to_csv('{}/{}_clinical_{}_life_table.csv'.format(outdir,basename,flag))
    return summary_score



def pre_data(project_dir):
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
    
    # ==== TCGA data
    tcga_AML_dir='{}/manuscript_figs/f11_clinical_TCGA_and_beatAML/TCGA_LAML/'.format(project_dir)
    tcga_AML_clinical_df = pd.read_csv('{}/f1_clinical_of_target_TF/selected_case_id_with_clinical_info.csv'.format(tcga_AML_dir),index_col=0) 
    tcga_AML_clinical_df = tcga_AML_clinical_df.fillna(0)
    tcga_AML_clinical_df['time_max']=tcga_AML_clinical_df[['days_to_last_follow_up', 'days_to_death']].max(axis=1)
    tcga_AML_clinical_df.loc[tcga_AML_clinical_df['vital_status']=='Dead','vital_status']=1
    tcga_AML_clinical_df.loc[tcga_AML_clinical_df['vital_status']=='Alive','vital_status']=0
    tcga_AML_expression_file='{}/f0_FPKM2TPM/TCGA_LAML_FPKM_Symbol.csv'.format(tcga_AML_dir)
    tcga_AML_expression_df = pd.read_csv(tcga_AML_expression_file,index_col=0)
    tcga_AML_expression_df = tcga_AML_expression_df[tcga_AML_clinical_df.index].dropna()

    return beatAML_expression_df, beatAML_clinical_df,tcga_AML_expression_df,tcga_AML_clinical_df



# == main ==

for gene_num in [9,31]:
    outdir='f5_cox_on_random_{}genes_beatAML_survival'.format(gene_num)
    os.makedirs(outdir,exist_ok=True)
    
    project_dir="/nv/vol190/zanglab/zw5j/wwork2018/AML_fran"
    project_dir="/Volumes/zanglab/zw5j/wwork2018/AML_fran"
    beatAML_expression_df,beatAML_clinical_df,tcga_AML_expression_df,tcga_AML_clinical_df = pre_data(project_dir)
    
    # ==== get the BART results for divergent genes
    # bart_results = '{}/updated_202010/f4_bart2_TR_expression/f2_PAML_DEG_DMG_removed_bart2/f1_DMG_removed_DEG_bart2_results/PAML_pthre5_rk3_ck2_genes_cluster3_div_bart_results.txt'.format(project_dir)
    # tfs_df = pd.read_csv(bart_results,sep='\t',index_col=0)
    # p_thre = 3
    # top_tfs = tfs_df[tfs_df['irwin_hall_pvalue']<np.power(.1,p_thre)].index;print(len(top_tfs))
    # print(len(top_tfs));continue
    
    # group ABC genes
    deg_clustering_dir = '{}/updated_202010/f1_kmeans_patients_by_deg/f1_kmeans_PAML_patients_by_PAML_DEG'.format(project_dir)
    basename = 'PAML_pthre5_rk3_ck2'
    deg_csv = deg_clustering_dir+'/csv/{}.csv'.format(basename)
    deg_list1 = [i.strip() for i in open(deg_clustering_dir+'/txt/{}_genes_cluster1.txt'.format(basename)).readlines()]
    deg_list2 = [i.strip() for i in open(deg_clustering_dir+'/txt/{}_genes_cluster2.txt'.format(basename)).readlines()]
    deg_list3 = [i.strip() for i in open(deg_clustering_dir+'/txt/{}_genes_cluster3.txt'.format(basename)).readlines()]
    deg_dict = {'groupA':deg_list1,
                'groupB':deg_list3,
                'groupC':deg_list2}
    
    
    # ==== repare the matrix for cox regression
    # penalizers = [.001, .002, .005, .01, .02, .05, .1, .2, .5]
    # l1_ratios = [0, .2, .4, .6, .8, 1]
    penalizers = [.08]
    l1_ratios = [1]
    # gene_num=9
    summary_score=pd.DataFrame()
    rand_gene_df = pd.DataFrame(index = np.arange(gene_num))
    
    for group in ['groupA','groupB','groupC'][:]:
        for rand_time in np.arange(5):
            group_genes = beatAML_expression_df.index.intersection(deg_dict[group])
            random_genes = np.random.choice(group_genes,gene_num,replace=False)
            beatAML_data=np.transpose(beatAML_expression_df.loc[random_genes].dropna()) # gene on the column
            cph_data = pd.concat([beatAML_clinical_df,beatAML_data],axis=1)
            # save the random selected genes
            rand_gene_df['{}_random{}'.format(group,rand_time)] = random_genes
            # use para of p = 0.001 and l1 = 0.4
            for penalizer in penalizers[:]:
                for l1_ratio in l1_ratios[:]:
                    basename = '{}_random{}_Cox_penalizer{}_l1_ratios{}'.format(group,rand_time,'p'.join(str(penalizer).split('.')),'p'.join(str(l1_ratio).split('.')))
            
                    ## ==== cross valication to test the parameters
                    k = 5;rerun_times = 10
                    cph = CoxPHFitter(penalizer=penalizer, l1_ratio=l1_ratio) # sparse solutions,
                    for ii in np.arange(rerun_times):
                        scores = k_fold_cross_validation(cph, cph_data,'overallSurvival',event_col='vitalStatus',k=k,scoring_method="concordance_index")
                        summary_score.loc[basename,'cv_time{}_mean'.format(ii)]=np.mean(scores).round(4)
                    cv_cols = [i for i in summary_score.columns if re.search('cv_time',i)]
                    summary_score.loc[basename,'cv_all_mean'] = summary_score.loc[basename,cv_cols].mean()
                    summary_score.loc[basename,'cv_all_max'] = summary_score.loc[basename,cv_cols].max()
            
                    # cph = CoxPHFitter()
                    cph = CoxPHFitter(penalizer=penalizer, l1_ratio=l1_ratio) # sparse solutions,
                    cph.fit(cph_data, duration_col='overallSurvival', event_col='vitalStatus')
                    cph.summary.to_csv(outdir+os.sep+'{}_cph_data.csv'.format(basename))
                    # print(cph.summary)
                    # ==== plot log HR
                    selected_genes = cph.summary['coef'][np.abs(cph.summary['coef'])>0].sort_values(ascending=False).index
                    plt.figure(figsize=(3,7))
                    # cph.plot()
                    cph.plot(selected_genes)
                    # plt.title('TR pthre {}'.format(np.power(.1,p_thre).round(3)),fontsize=16)
                    plt.ylim([-1,len(selected_genes)])
                    plt.title(basename.split('Cox_')[1])
                    plt.savefig(outdir+os.sep+'{}_cph_logHR.pdf'.format(basename),bbox_inches='tight',pad_inches=.1,dpi=600,transparent=True)
                    plt.show()
                    plt.close()
                    # print(len(selected_genes))
            
                    ## ==== survival anlaysis of beatAML data
                    risk_days = [0,250,500,750,1000,1250,1500]
                    logrank_clinical_cols = ['vitalStatus','overallSurvival']
                    calinical_analysis(summary_score,cph,beatAML_data,beatAML_clinical_df,logrank_clinical_cols,'BeatAML',basename,risk_days)
                        
    ## == save out the results
    summary_score.to_csv('{}/summary_score.csv'.format(outdir))        
    rand_gene_df.to_csv('{}/random_genes.csv'.format(outdir))      

