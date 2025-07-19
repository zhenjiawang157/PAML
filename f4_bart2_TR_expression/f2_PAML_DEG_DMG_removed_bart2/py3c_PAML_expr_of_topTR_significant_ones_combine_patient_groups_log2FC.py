import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
#from GenomeData import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=15
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
import seaborn as sns
# sns.set(font_scale=1.2)
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
    log2FC = pd.read_csv(file_dir+os.sep+'{}_deseq2_log2FC.txt'.format(cohort),sep=',',index_col=0)
    
    return log2FC



def expr_plot_for_patients(top_tfs,patients1,patients2,log2FC,cohort,group,outdir):
    
    # ==== only plot those TFs with sig diff expression
    tf_expr_compr_ttest1 = pd.read_csv('f3a_PAML_expr_of_topTR_sep_patient'+os.sep+'ttest_{}_{}_patientG1.csv'.format(cohort,group),index_col=0)
    tf_expr_compr_ttest2 = pd.read_csv('f3a_PAML_expr_of_topTR_sep_patient'+os.sep+'ttest_{}_{}_patientG2.csv'.format(cohort,group),index_col=0)
    g1_sig_tfs = tf_expr_compr_ttest1[tf_expr_compr_ttest1['pvalue']<0.05].index
    g2_sig_tfs = tf_expr_compr_ttest2[tf_expr_compr_ttest2['pvalue']<0.05].index
    top_tfs = [i for i in top_tfs if i in g1_sig_tfs.union(g2_sig_tfs).union(['KMT2A'])][::-1]
    print(top_tfs)
    
    box_vals1,box_vals2,positions,gene_names=[],[],np.array([]),[]
    face_colors1,face_colors2=[],[]
    df_patient_c1,df_patient_c2 = pd.DataFrame(),pd.DataFrame()
    for ii in np.arange(len(top_tfs)):
        positions = np.append(positions,ii)
        # alias of some genes
        gene_name = top_tfs[ii]
        if gene_name=='H2AZ':
            gene_name='H2AFZ'
        if gene_name=='KMT2A':
            gene_name='MLL'
         
        face_colors1.append('tab:orange' if gene_name in g1_sig_tfs else 'w')
        face_colors2.append('tab:blue' if gene_name in g2_sig_tfs else 'w')
        
        box_vals1.append(log2FC.loc[gene_name][patients1].values)
        box_vals2.append(log2FC.loc[gene_name][patients2].values)
        gene_names.append(gene_name)

        # to save the log2FC values
        df_patient_c1 = pd.concat([df_patient_c1,log2FC.loc[gene_name][patients1].to_frame().rename(columns={0:gene_name})],axis=1)
        df_patient_c2 = pd.concat([df_patient_c2,log2FC.loc[gene_name][patients2].to_frame().rename(columns={0:gene_name})],axis=1)
    
    # use AML_xxx format as index
    df_patient_c1.index = info.loc[patients1]['subject_ID']
    df_patient_c2.index = info.loc[patients2]['subject_ID']
    # reverse the columns
    df_patient_c1 = df_patient_c1[df_patient_c1.columns[::-1]]
    df_patient_c2 = df_patient_c2[df_patient_c2.columns[::-1]]
    # save the df
    df_patient_c1.to_csv('{}/{}_genes_{}_patient_C1.csv'.format(outdir,cohort,group))
    df_patient_c2.to_csv('{}/{}_genes_{}_patient_C2.csv'.format(outdir,cohort,group))
    
    # plot the figs
    face_colors1 = ['tab:orange']*top_tf_num
    face_colors2 = ['tab:blue']*top_tf_num
    
    figname = '{}/{}_genes_{}.pdf'.format(outdir,cohort,group)
    plt.figure(figsize=(6,3))
    # boxes of patient group1
    g = plt.boxplot(box_vals1,positions=positions+.15,vert=False,widths = .2,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=True,lw=1),\
                medianprops=dict(color='grey'),showfliers=False)      
    for patch, color in zip(g['boxes'], face_colors1):
        patch.set_facecolor(color)
        # patch.set_edgecolor(color)

    # boxes of patient group2
    g = plt.boxplot(box_vals2,positions=positions-.15,vert=False,widths = .2,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=True,lw=1),\
                medianprops=dict(color='grey'),showfliers=False)      
    for patch, color in zip(g['boxes'], face_colors2):
        patch.set_facecolor(color)


    if group=='div':
        fc_thre=6
        barh_pos=[1,3.,5.]
    if group=='down':
        fc_thre=5
        barh_pos=[0.,2.,4.]
    if group=='up':
        fc_thre=2
        barh_pos=[1.,3.]
        
    plt.barh(barh_pos,fc_thre,height=1,color=['gainsboro']*3,edgecolor = "w",lw=0,zorder=0,clip_on=False,alpha=.7)
    plt.barh(barh_pos,-1*fc_thre,height=1,color=['gainsboro']*3,edgecolor = "w",lw=0,zorder=0,clip_on=False,alpha=.7)
    
    plt.axes().set_yticks(np.arange(len(top_tfs)))
    plt.axes().set_yticklabels(top_tfs,rotation=0,ha='right')
    plt.axes().set_xticks(np.arange(-1*fc_thre,fc_thre+1))
#     plt.legend([g["boxes"][0],g["boxes"][1]],['Diagnosis','Relapse'],fontsize=12,borderaxespad=0.2,labelspacing=.2,handletextpad=0.2,handlelength=1,loc=0,frameon=False)
    plt.xlabel('log2 foldchange',)
    plt.ylim([-.5,len(top_tfs)-.5])
    plt.xlim([-1*fc_thre,fc_thre])
    plt.axes().axvline(x=0,color='k')
#     plt.title('{} {}'.format(cohort,group))
    plt.show()
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.close()




# ==== main section
outdir = 'f3c_PAML_expr_of_topTR_combine_patient_sig_log2FC'
os.makedirs(outdir,exist_ok=True)

project_dir="/nv/vol190/zanglab/zw5j/wwork2018/AML_fran"
# project_dir="/Volumes/zanglab/zw5j/wwork2018/AML_fran"
patient_clustering_dir = '{}/updated_202010/f2_kmeans_patients_by_div_gene/f1_kmeans_PAML_patients_by_PAML_div_genes'.format(project_dir)
clinical_dir = '{}/f0_fran_data_new/data_processed/'.format(project_dir)

for group in ['div','up','down'][:]:
    for cohort in ["PAML"]:
        basename = '{}_{}'.format(cohort,group)
        info = pd.read_csv(clinical_dir+os.sep+'{}_clinical.csv'.format(cohort),index_col=0)

        #==== bart top tf
        top_tf_num = 10
        bart_results_file=glob.glob('f1_DMG_removed_DEG_bart2_results/PAML_pthre5_rk3_ck2_genes_*{}_bart_results.txt'.format(group))[0]
        bart_results = pd.read_csv(bart_results_file,sep='\t',index_col=0)
        top_tfs = bart_results.index[:top_tf_num]
        

        # ==== gene expression TPM
        # tpm_diagnosis,tpm_relapse =  return_patient_info(project_dir,cohort)
        log2FC =  return_patient_info(project_dir,cohort)
        
        basename = '{}_pthre5_rk3_ck2'.format(cohort)
        patient_list1_tmp = [i.strip() for i in open(patient_clustering_dir+'/txt/{}_patients_cluster1.txt'.format(basename)).readlines()]
        patient_list2_tmp = [i.strip() for i in open(patient_clustering_dir+'/txt/{}_patients_cluster2.txt'.format(basename)).readlines()]

        # ==== re-order the patients followed by down/up-genes
        # re-rank the patient clusters
        patient_list1,patient_list2 = patient_list2_tmp,patient_list1_tmp
        expr_plot_for_patients(top_tfs,patient_list1,patient_list2,log2FC,cohort,group,outdir)

           
        
        
        
 




