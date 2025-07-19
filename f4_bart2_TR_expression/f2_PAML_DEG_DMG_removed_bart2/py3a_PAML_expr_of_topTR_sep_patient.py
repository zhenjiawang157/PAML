import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
#from GenomeData import *
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")
from scipy import stats
#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')


def return_patient_info(project_dir,cohort):
    
    '''This is used to return the deseq2 results and TPM for patient samples '''
    file_dir = '{}/f0_fran_data_new/data_processed'.format(project_dir)
    
    tpm_diagnosis = pd.read_csv(file_dir+os.sep+'{}_TPM_diagnosis.txt'.format(cohort),sep=',',index_col=0)
    tpm_relapse = pd.read_csv(file_dir+os.sep+'{}_TPM_relapse.txt'.format(cohort),sep=',',index_col=0)
    tpm_normal = pd.read_csv(file_dir+os.sep+'{}_TPM_normal.txt'.format(cohort),sep=',',index_col=0)
    
    # return tpm_diagnosis,tpm_relapse
    return np.log2(tpm_diagnosis + 1),np.log2(tpm_relapse + 1)

    
    

def mark_pvalue(compr_pos,positions,box_vals):
    s,p = stats.ttest_rel(box_vals[compr_pos[1]],box_vals[compr_pos[0]],nan_policy='omit')
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),96)*1.00 ,1.05, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),0)*0.99
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    p_label='*'
    if p<0.05:
#         if p<0.001:
#             p_label = '**'
        if compr_pos[2] == 't':
#             plt.plot([x1*1.01, x1*1.01, x2*.99, x2*.99], [y, y*h, y*h, y], lw=1, c=col)
            plt.text((x1+x2)*.5, y*.95, p_label, ha='center', va='bottom', color=col,fontsize=16)
        else:
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y2, y2-2, y2-2, y2], lw=1, c=col)
            plt.text((x1+x2)*.5, y2-2, p_label, ha='center', va='top', color=col,fontsize=16)
    return s,p


def compr_boxplot(positions,box_vals,figname,top_tfs,colors,cohort,group,flag):
    
    plt.figure(figsize=(5,3.2))
    g = plt.boxplot(box_vals,positions=positions,widths = .15,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=True,lw=1),\
                medianprops=dict(color='grey'),showfliers=False)    
    for patch, color in zip(g['boxes'], colors):
        patch.set_facecolor(color)

#     scatter_X = []
#     for position_id in np.arange(len(positions)):
#         scatter_x = np.random.normal(positions[position_id],0.07,len(box_vals[position_id]))
#         plt.scatter(scatter_x,box_vals[position_id],color=colors[position_id],s=10,zorder=0,alpha=0.99)
    
    tf_expr_compr_ttest = pd.DataFrame()
    for i in np.arange(0,len(top_tfs)):
        for compr_pos in [[2*int(i),2*int(i)+1,'t']]:
            s,p = mark_pvalue(compr_pos,positions,box_vals)
            tf_expr_compr_ttest.loc[top_tfs[i],'t-statistic'] = s
            tf_expr_compr_ttest.loc[top_tfs[i],'pvalue'] = p

    plt.axes().set_xticks(np.arange(len(top_tfs))+.17)
    plt.axes().set_xticklabels(top_tfs,fontsize=16,rotation=45,ha='right')
    plt.legend([g["boxes"][0],g["boxes"][1]],['Diagnosis','Relapse'],fontsize=12,borderaxespad=0.2,labelspacing=.2,handletextpad=0.2,handlelength=1,loc=0,frameon=False)
    plt.ylabel('log2 TPM+1',fontsize=16)
    plt.xlim([-.5,len(top_tfs)-.2])
    plt.title('{} {} {}'.format(cohort,group,flag))
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.show()
    plt.close()
    return tf_expr_compr_ttest


def expr_plot_for_patients(top_tfs,patients,tpm_diagnosis,tpm_relapse,cohort,group,flag,outdir):
    
    box_vals,positions,gene_names=[],[],[]
    for ii in np.arange(len(top_tfs)):
        positions.append(ii)
        # alias of some genes
        gene_name = top_tfs[ii]
        if gene_name=='H2AZ':
            gene_name='H2AFZ'
        if gene_name=='KMT2A':
            gene_name='MLL'
            
        box_vals.append(tpm_diagnosis.loc[gene_name][patients].values)
        # relapse
        positions.append(ii+.25)
        box_vals.append(tpm_relapse.loc[gene_name][patients].values)
        gene_names.append(gene_name)
    colors = ['b','r']*top_tf_num
    figname = '{}/{}_{}_patient{}.pdf'.format(outdir,cohort,group,flag)
    tf_expr_compr_ttest = compr_boxplot(positions,box_vals,figname,gene_names,colors,cohort,group,flag)
    tf_expr_compr_ttest.to_csv(outdir+os.sep+'ttest_{}_{}_patient{}.csv'.format(cohort,group,flag))
    print('cohort {}, #{} {}'.format(cohort,flag,len(patient_list1)))



# ==== main section
outdir = 'f3a_PAML_expr_of_topTR_sep_patient_ttest_rel'
os.makedirs(outdir,exist_ok=True)

project_dir="/nv/vol190/zanglab/zw5j/wwork2018/AML_fran"
project_dir="/Volumes/zanglab/zw5j/wwork2018/AML_fran"
patient_clustering_dir = '{}/updated_202010/f2_kmeans_patients_by_div_gene/f1_kmeans_PAML_patients_by_PAML_div_genes'.format(project_dir)


for group in ['div','up','down'][:]:
    for cohort in ["PAML"]:
        basename = '{}_{}'.format(cohort,group)
        
        #==== bart top tf
        top_tf_num = 10
        bart_results_file=glob.glob('f1_DMG_removed_DEG_bart2_results/PAML_pthre5_rk3_ck2_genes_*{}_bart_results.txt'.format(group))[0]
        bart_results = pd.read_csv(bart_results_file,sep='\t',index_col=0)
        top_tfs = bart_results.index[:top_tf_num]
        
        # ==== gene expression TPM
        tpm_diagnosis,tpm_relapse =  return_patient_info(project_dir,cohort)
        print('\n',top_tfs.difference(tpm_diagnosis.index))
        # top_tfs = top_tfs.intersection(tpm_diagnosis.index)
        
        basename = '{}_pthre5_rk3_ck2'.format(cohort)
        patient_list1_tmp = [i.strip() for i in open(patient_clustering_dir+'/txt/{}_patients_cluster1.txt'.format(basename)).readlines()]
        patient_list2_tmp = [i.strip() for i in open(patient_clustering_dir+'/txt/{}_patients_cluster2.txt'.format(basename)).readlines()]

        # ==== re-order the patients followed by down/up-genes
        # re-rank the patient clusters
        patient_list1,patient_list2 = patient_list2_tmp,patient_list1_tmp

        # for patient group 1
        expr_plot_for_patients(top_tfs,patient_list1,tpm_diagnosis,tpm_relapse,cohort,group,'G1',outdir)

        # for patient group 1
        expr_plot_for_patients(top_tfs,patient_list2,tpm_diagnosis,tpm_relapse,cohort,group,'G2',outdir)
           
        
        
        
# gene = 'TAL1'
# a = tpm_diagnosis.loc[gene,patient_list2]
# b = tpm_relapse.loc[gene,patient_list2]
# stats.ttest_rel(a,b)

# x = log2FC.loc[gene,patient_list1]
# y = log2FC.loc[gene,patient_list2]
# pvalue = pd.read_csv(file_dir+os.sep+'{}_deseq2_qvalue.txt'.format(cohort),sep=',',index_col=0)
# pvalue.loc[gene,patient_list2]

# plt.boxplot([x,y])
# plt.boxplot([np.log2(a),np.log2(b)])
