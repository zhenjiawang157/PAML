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
    plt.title('{} {}'.format(basename,flag),fontsize=15)
    plt.xlabel('Days',fontsize=15)
    plt.ylabel('Survival probability',fontsize=15)
    plt.xticks(risk_days)
    plt.legend(bbox_to_anchor=[1.0,1.0],loc='upper right',frameon=False,\
               borderaxespad=0.2,labelspacing=.2,handletextpad=0.2,handlelength=1)
    plt.savefig(figname,bbox_inches='tight',pad_inches=.1,dpi=300,transparent=True)
    plt.show()
    plt.close()
    
    return results







# == main ==
outdir='f3_survival_TR9'
os.makedirs(outdir,exist_ok=True)

# project_dir="/nv/vol190/zanglab/zw5j/wwork2018/AML_fran"
project_dir="/Volumes/zanglab/zw5j/wwork2018/AML_fran"
beatAML_dir='{}/manuscript_figs/f11_clinical_TCGA_and_beatAML/beatAML/'.format(project_dir)
beatAML_RPKM = pd.read_csv('{}/beatAML_RPKM.csv'.format(beatAML_dir),index_col=0)

beatAML_clinical_df = pd.read_excel('{}/beatAML_suppl_ClinicalInformation_7+3Tx.xlsx'.format(beatAML_dir),index_col=0,sheet_name='7+3') 
beatAML_clinical_df = beatAML_clinical_df[['overallSurvival','vitalStatus']].dropna()
# ==== check the patients with both clinical data and expression data
shared_patients = beatAML_clinical_df.index.intersection(beatAML_RPKM.columns)
beatAML_RPKM = beatAML_RPKM[shared_patients]
beatAML_clinical_df = beatAML_clinical_df.loc[shared_patients]
beatAML_clinical_df.loc[beatAML_clinical_df['vitalStatus']=='Dead','vitalStatus']=1
beatAML_clinical_df.loc[beatAML_clinical_df['vitalStatus']=='Alive','vitalStatus']=0
beatAML_clinical_df.columns = ['time','status']
logrank_clinical_data = beatAML_clinical_df
# logrank_clinical_data.to_csv('{}/clinical_data.csv'.format(outdir))      


# selected_TRs = ['GATA1','RAD21','CUX1','KLF16','ADNP','ZEB2','TEAD4','STAG1','SIRT6']
coef = pd.read_csv('f2_figs/coef_s0p08.csv',header=None,index_col=0)
coef = coef[coef!=0].dropna()
coef.columns = ['coef']
coef.to_csv('{}/coef.csv'.format(outdir))
expr_data = np.transpose(beatAML_RPKM.loc[coef.index])
inner_product = np.inner(expr_data,coef['coef'])
inner_df = pd.DataFrame(inner_product)
inner_df.index = expr_data.index
c1_lower = inner_df.sort_values(by=[0],ascending=True).index[:int(.5*inner_df.shape[0])]
c2_higher = inner_df.sort_values(by=[0],ascending=True).index[-1*int(.5*inner_df.shape[0]):]

# plot the clinical figs
basename = 'TR9'
flag = 'beatAML'
risk_days = [0,250,500,750,1000,1250,1500]
summary_score=pd.DataFrame()

survival_df = pd.DataFrame(index = np.append(c1_lower,c2_higher))#;print(survival_df)
survival_df.loc[set(c1_lower).intersection(logrank_clinical_data.index),'group']='c1'
survival_df.loc[set(c2_higher).intersection(logrank_clinical_data.index),'group']='c2'
survival_df['status']= logrank_clinical_data.loc[survival_df.index]['status']         
survival_df['time']= logrank_clinical_data.loc[survival_df.index]['time']
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

summary_score.to_csv('{}/summary_score.csv'.format(outdir))        


