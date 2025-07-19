import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
#from GenomeData import *
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=9
import seaborn as sns
sns.set(font_scale=1)
sns.set_style("whitegrid", {'axes.grid' : False})
from matplotlib import gridspec
#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
from  scipy.cluster.hierarchy import fcluster

def expr_with_info_plot(df,counts_row,counts_col,mutation,cytogenetic,outdir,basename,metric='euclidean',method='average'):

    
    plt.figure(figsize = (4,8))
    height_ratios = [5,1.3,1.7]
    # width_ratios=[5,.3]
    gs = gridspec.GridSpec(3,1,height_ratios=height_ratios,hspace=0.04) 
    

    ## ==== heatmap of gene expression data
    loc=0  
    ax = plt.subplot(gs[loc,0])       
    g = sns.heatmap(df,cmap=plt.cm.PiYG_r,cbar=False,vmax=2,vmin=-2\
                    ,xticklabels=False,yticklabels=False,cbar_kws={"orientation": "vertical","use_gridspec":False}\
                    ,ax=ax,rasterized=True)#"use_gridspec":False,"location":"top"})    
#     ax.xaxis.set_ticks_position('top')
#     ax.set_xticklabels(g.get_xticklabels(),size = 8,rotation=90)
    ax.set_xlabel(None)
    ax.set_ylabel('\n\n{} genes'.format(df.shape[0]),fontsize=14)
    ax.yaxis.set_label_position('left')
    
    for i in np.arange(len(counts_row)+1):
        ax.axhline(y=sum(counts_row[:i]),xmin=-.00,xmax=1.00,color='k',lw=1.5,clip_on=False)
    for i in np.arange(len(counts_col)+1):
        ax.axvline(x=sum(counts_col[:i]),ymin=-.00,ymax=1.00,color='k',lw=1.5,clip_on=False)
    

    ## ==== mutation matrix, keep those genes w/ > mutation in current patient
    extra_patients = df.columns.difference(mutation.index)
    mutation = pd.concat([mutation,pd.DataFrame(index = extra_patients)])
    df_mutation = mutation.loc[df.columns]
    df_mutation = np.transpose(df_mutation)
    df_mutation = df_mutation.loc[selected_genes];
    print(selected_genes,df_mutation)
    
    loc=1
    ax = plt.subplot(gs[loc,0]) 
    g1 = sns.heatmap(df_mutation,cmap=plt.cm.Greys,xticklabels=False, yticklabels=True,\
                    linewidths=.05,linecolor='k',\
                    cbar=False,vmin=0,vmax=1,ax = ax)#,cbar_kws={"orientation": "horizontal"})#"use_gridspec":False,"location":"top"})    
    ax.set_xticklabels(g1.get_xticklabels(),size = 8)
    ax.set_yticklabels(g1.get_yticklabels(),size = 8)
    ax.set_xlabel(None)


    ## ==== cytogenetic matrix, keep those genes w/ > mutation in current patient
    extra_patients = df.columns.difference(cytogenetic.index)
    cytogenetic = pd.concat([cytogenetic,pd.DataFrame(index = extra_patients)])
    cytogenetic = cytogenetic.loc[df.columns]
    cytogenetic = np.transpose(cytogenetic)

    loc=2
    ax = plt.subplot(gs[loc,0]) 
    g1 = sns.heatmap(cytogenetic,cmap=plt.cm.Greys,xticklabels=True, yticklabels=True,\
                    linewidths=.05,linecolor='k',\
                    cbar=False,vmin=0,vmax=1,ax = ax)#,cbar_kws={"orientation": "horizontal"})#"use_gridspec":False,"location":"top"})    
    ax.set_xticklabels(g1.get_xticklabels(),size = 8)
    ax.set_yticklabels(g1.get_yticklabels(),size = 8)
    ax.set_xlabel(None)
    
    plt.savefig(outdir+'/{}.pdf'.format(basename),bbox_inches = 'tight',pad_inches=0.1,transparent=True,dpi=600)
    plt.show()
    plt.close()





## ==== main section
    
outdir = 'f2_SG_heatmap_with_mutation'
os.makedirs(outdir,exist_ok=True)

project_dir="/nv/vol190/zanglab/zw5j/wwork2018/AML_fran"
project_dir="/Volumes/zanglab/zw5j/wwork2018/AML_fran"


# ==== mutation info, genes with # mutation >10 across all PAML patients, not limited to 59 paired
mutation_df = pd.read_csv('{}/ff5_CNV_mutation/f2_mutation/PAML_gene_mutation_relapse_GT10.csv'.format(project_dir),index_col=0)
# mutation=1.1, non-mutation=0.1, not-detected = NA
mutation_df = mutation_df.fillna(0)+0.2 
selected_genes = ['DNMT3A','FLT3_ITD','NPM1','NRAS','BRAF',\
                  'SETBP1','TET2','WT1','ASXL1','CEBPA']
mutation_df = mutation_df[selected_genes]

# ==== cytogenetics 
cytogenetic_df = pd.read_excel('{}/ff5_CNV_mutation/data/PAML_cytogenetics_category_col_F-S.xlsx'.format(project_dir),index_col=0)
# selected_cols=['NCG','Other','Complex']
selected_cols=['NCG','Complex','Inv(16)','8+','t(8;21)','del5(q)','7-','17-/17p abn','11q23','5-','t(9;22)','t(6;9)','Inv(3) or t(3;3)']
cytogenetic_df = cytogenetic_df.fillna(0)+0.2 
cytogenetic_df = cytogenetic_df[selected_cols]


# dir of DEG and clinical info
# use PAML DEG
deg_clustering_dir = '{}/updated_202010/f1_kmeans_patients_by_deg/f1_kmeans_PAML_patients_by_PAML_DEG'.format(project_dir)
patient_clustering_dir = '{}/updated_202010/f2_kmeans_patients_by_div_gene/f3_kmeans_SG_patients_by_PAML_div_SG_DEG'.format(project_dir)
clinical_dir = '{}/f0_fran_data_new/data_processed/'.format(project_dir)

for cohort in ['SG']:
    
    # ==== clinical match ID
    info = pd.read_csv(clinical_dir+os.sep+'{}_clinical.csv'.format(cohort),index_col=0)
    # ==== clinical time-to-relapse analysis for each gene
    # clinical_df = pd.read_csv('{}/ff2_survival/f5_survival_each_gene_by_TPM/{}_survival_results.csv'.format(project_dir,cohort),index_col=0)
    
    # ==== info of deg clustering
    basename = '{}_pthre5_rk3_ck2'.format(cohort)
    deg_csv = patient_clustering_dir+'/csv/{}.csv'.format(basename)
    deg_df = pd.read_csv(deg_csv,index_col=0)
    # PAML DEG
    deg_list1 = [i.strip() for i in open(deg_clustering_dir+'/txt/PAML_pthre5_rk3_ck2_genes_cluster1.txt'.format(basename)).readlines()]
    deg_list2 = [i.strip() for i in open(deg_clustering_dir+'/txt/PAML_pthre5_rk3_ck2_genes_cluster2.txt'.format(basename)).readlines()]
    deg_list3 = [i.strip() for i in open(deg_clustering_dir+'/txt/PAML_pthre5_rk3_ck2_genes_cluster3.txt'.format(basename)).readlines()]
    # keep only those SG DEG
    deg_list1 = [i for i in deg_list1 if i in deg_df.index]
    deg_list2 = [i for i in deg_list2 if i in deg_df.index]
    deg_list3 = [i for i in deg_list3 if i in deg_df.index]
    with open(outdir+os.sep+'{}_gene1.txt'.format(basename),'w') as outf:
        outf.write('\n'.join(deg_list1)+'\n')
    with open(outdir+os.sep+'{}_gene2.txt'.format(basename),'w') as outf:
        outf.write('\n'.join(deg_list2)+'\n')
    with open(outdir+os.sep+'{}_gene3.txt'.format(basename),'w') as outf:
        outf.write('\n'.join(deg_list3)+'\n')
    
    patient_list1 = [i.strip() for i in open(patient_clustering_dir+'/txt/{}_patients_cluster1.txt'.format(basename)).readlines()]
    patient_list2 = [i.strip() for i in open(patient_clustering_dir+'/txt/{}_patients_cluster2.txt'.format(basename)).readlines()]
    patient_list1_new = info.loc[patient_list1]['subject_ID']
    patient_list2_new = info.loc[patient_list2]['subject_ID']
    with open(outdir+os.sep+'{}_patient1.txt'.format(basename),'w') as outf:
        outf.write('\n'.join(patient_list1)+'\n')
    with open(outdir+os.sep+'{}_patient2.txt'.format(basename),'w') as outf:
        outf.write('\n'.join(patient_list2)+'\n')
    with open(outdir+os.sep+'{}_patient1_subjectID.txt'.format(basename),'w') as outf:
        outf.write('\n'.join(patient_list1_new)+'\n')
    with open(outdir+os.sep+'{}_patient2_subjectID.txt'.format(basename),'w') as outf:
        outf.write('\n'.join(patient_list2_new)+'\n')
        
    # ==== re-order the genes followed by down-/div-/up-genes
    reindex = deg_list1+deg_list3+deg_list2
    recol = patient_list2+patient_list1

    # position for cluster box
    counts_row = [len(deg_list1),len(deg_list3),len(deg_list2)]
    counts_col = [len(patient_list2),len(patient_list1)]
        
    deg_df = np.transpose(np.transpose(deg_df.reindex(reindex)).reindex(recol))

    # ==== convert the ID into AML_xxx format
    # info_df = info.loc[deg_df.columns][['time_until_relapse_days']].fillna(0);print(info_df)
    # info_df = ((info_df/365).astype(int)+1).clip(upper=3)
    deg_df.columns = info.loc[deg_df.columns]['subject_ID']
    
#     for metric in ['euclidean','cosine','correlation']:
#         for method in ['average','single','weighted','ward','complete','centroid']:
    expr_with_info_plot(deg_df,counts_row,counts_col,mutation_df,cytogenetic_df,outdir,basename)
    

    # calculate the % down-groupA genes and up-groupC genes
    df_sub = deg_df.loc[deg_list1]
    percentage_each_gene = 100*(df_sub<0).sum().sort_values(ascending=False)/df_sub.shape[0]
    pd.DataFrame(percentage_each_gene,columns=['% down genes']).round(2).to_csv(outdir+os.sep+'{}_percent_groupA_down_each_patient.csv'.format(basename))
    percentage_each_patient = 100*(df_sub<0).sum(axis=1).sort_values(ascending=False)/df_sub.shape[1]
    pd.DataFrame(percentage_each_patient,columns=['% down patients']).round(2).to_csv(outdir+os.sep+'{}_percent_groupA_down_each_gene.csv'.format(basename))
    print((percentage_each_patient>50).sum()/percentage_each_patient.shape[0])
    
    df_sub = deg_df.loc[deg_list2]
    percentage_each_gene = 100*(df_sub>0).sum().sort_values(ascending=False)/df_sub.shape[0]
    pd.DataFrame(percentage_each_gene,columns=['% up genes']).round(2).to_csv(outdir+os.sep+'{}_percent_groupC_up_each_patient.csv'.format(basename))
    percentage_each_patient = 100*(df_sub>0).sum(axis=1).sort_values(ascending=False)/df_sub.shape[1]
    pd.DataFrame(percentage_each_patient,columns=['% up patients']).round(2).to_csv(outdir+os.sep+'{}_percent_groupC_up_each_gene.csv'.format(basename))
    print((percentage_each_patient>50).sum()/percentage_each_patient.shape[0])