import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
#from GenomeData import *
import matplotlib
matplotlib.use('Agg')
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

def expr_with_info_plot(df,counts_row,col_colors,mutation,cytogenetic,outdir,basename,metric='euclidean',method='average'):

    
    plt.figure(figsize = (4,6))
    height_ratios = [5,1.5,1.5]
    gs = gridspec.GridSpec(3,1,height_ratios=height_ratios,hspace=0.02) 
    
    ## ==== heatmap of gene expression data
#    loc=0  
#    ax = plt.subplot(gs[loc,0])       
#    g = sns.heatmap(df,cmap=plt.cm.PiYG_r,cbar=False,vmax=2.5,vmin=-2.5\
#                    ,xticklabels=True,yticklabels=False,cbar_kws={"orientation": "vertical","use_gridspec":False})#"use_gridspec":False,"location":"top"})    
#    ax.xaxis.set_ticks_position('top')
#    ax.set_xticklabels(g.get_xticklabels(),size = 8,rotation=90)
#    ax.set_xlabel(None)
#    ax.set_ylabel('{} Genes'.format(df.shape[0]),fontsize=16)

     # ==== if you need to plot using clustermap
    g = sns.clustermap(df,figsize=(4,8),cmap=plt.cm.PiYG_r,vmax=2,vmin=-2\
                     ,xticklabels=False,yticklabels=False,row_cluster=False\
                     ,metric=metric,method=method\
                     ,col_colors = col_colors\
                     ,cbar_kws={"orientation": "vertical","use_gridspec":False},rasterized=True)
    # g.cax.set_position([.15,.0,0.03,.2])
    g.cax.set_ylabel('log$_2$ FC')
    g.gs.update(top=0.95, bottom=0.41,left=-0.075)
    g.ax_heatmap.set_xlabel(None)
    g.ax_heatmap.set_ylabel('{} genes'.format(df.shape[0]),fontsize=12)
    g.cax.set_visible(False)
    g.ax_heatmap.set_title('metric = {}\nmethod = {}\n\n\n\n'.format(metric,method),fontsize=18)
    # for i in np.arange(len(counts_row)+1):
    #     plt.axes().axhline(y=sum(counts_row[:i]),xmin=-.00,xmax=1.00,color='k',lw=1.5,clip_on=False)
    
    
    ## === save out the patient samples
    clustered_patients = df.columns[g.dendrogram_col.dendrogram['leaves']]
    z = g.dendrogram_col.linkage # cluster assign
    cluster_num = 2
    cluster_def = fcluster(z,t=cluster_num,criterion='maxclust')
    for i in np.arange(cluster_num):
        cluster_ids = df.columns[cluster_def==i+1]
        with open(outdir+'/txt/{}_{}_{}_cluster{}.txt'.format(basename,metric,method,i+1),'w') as outf:
            outf.write('\n'.join(cluster_ids)+'\n')

    #### mutation matrix, keep those genes w/ > mutation in current patient
    extra_patients = df.columns.difference(mutation.index)
    mutation = pd.concat([mutation,pd.DataFrame(index = extra_patients)])
    df_mutation = mutation.loc[clustered_patients]
    df_mutation = np.transpose(df_mutation)
    df_mutation = df_mutation.loc[selected_genes];
    
    loc=1
    ax = g.fig.add_subplot(gs[loc,0])    
    g1 = sns.heatmap(df_mutation,cmap=plt.cm.Greys,xticklabels=False, yticklabels=True,\
                    linewidths=.05,linecolor='k',\
                    cbar=False,vmin=0,vmax=1,ax = ax)#,cbar_kws={"orientation": "horizontal"})#"use_gridspec":False,"location":"top"})    
    ax.set_xticklabels(g1.get_xticklabels(),size = 8)
    ax.set_yticklabels(g1.get_yticklabels(),size = 8)
    ax.set_xlabel(None)


    #### cytogenetic matrix, keep those genes w/ > mutation in current patient
    extra_patients = df.columns.difference(cytogenetic.index)
    cytogenetic = pd.concat([cytogenetic,pd.DataFrame(index = extra_patients)])
    cytogenetic = cytogenetic.loc[clustered_patients]
    cytogenetic = np.transpose(cytogenetic)
    cytogenetic = cytogenetic.loc[(cytogenetic>.2).sum(axis=1)>0]

    loc=2
    ax = g.fig.add_subplot(gs[loc,0])    
    g1 = sns.heatmap(cytogenetic,cmap=plt.cm.Greys,xticklabels=True, yticklabels=True,\
                    linewidths=.05,linecolor='k',\
                    cbar=False,vmin=0,vmax=1,ax = ax)#,cbar_kws={"orientation": "horizontal"})#"use_gridspec":False,"location":"top"})    
    ax.set_xticklabels(g1.get_xticklabels(),size = 8)
    ax.set_yticklabels(g1.get_yticklabels(),size = 8)
    ax.set_xlabel(None)

    
    plt.show()
    plt.savefig(outdir+'/{}_{}_{}.pdf'.format(basename,metric,method),bbox_inches = 'tight',pad_inches=0.1,transparent=True,dpi=600)
    plt.close()





### main section
    
outdir = 'f5_PAML_HC_heatmap_with_mutation'
os.makedirs(outdir,exist_ok=True)
os.makedirs(outdir+os.sep+'txt',exist_ok=True)

project_dir="/nv/vol190/zanglab/zw5j/wwork2018/AML_fran"
# project_dir="/Volumes/zanglab/zw5j/wwork2018/AML_fran"


# ==== mutation info, genes with # mutation >10 across all PAML patients, not limited to 59 paired
mutation_df = pd.read_csv('{}/ff5_CNV_mutation/f2_mutation/PAML_gene_mutation_relapse_GT10.csv'.format(project_dir),index_col=0)
# mutation=1.1, non-mutation=0.1, not-detected = NA
mutation_df = mutation_df.fillna(0)+0.2 
selected_genes = ['DNMT3A','FLT3_ITD','NPM1','NRAS','BRAF',\
                  'SETBP1','TET2','WT1','ASXL1','CEBPA']
mutation_df = mutation_df[selected_genes]

# ==== cytogenetics 
cytogenetic_df = pd.read_excel('{}/ff5_CNV_mutation/data/PAML_cytogenetics_category_col_F-S.xlsx'.format(project_dir),index_col=0)
selected_cols=['NCG','Other','Complex']
selected_cols=['NCG','Complex','Inv(16)','8+','t(8;21)','del5(q)','7-','17-/17p abn','11q23','5-','t(9;22)','t(6;9)','Inv(3) or t(3;3)','Other']
cytogenetic_df = cytogenetic_df.fillna(0)+0.2 
cytogenetic_df = cytogenetic_df[selected_cols]


# dir of DEG and clinical info
# deg_clustering_dir = '{}/manuscript_figs/f1_DEG_with_mutation/f4_kmeans_clustering_by_DIV___DEL'.format(project_dir)
deg_clustering_dir = '{}/updated_202010/f1_kmeans_patients_by_deg/f1_kmeans_PAML_patients_by_PAML_DEG'.format(project_dir)
patient_clustering_dir = '{}/updated_202010/f2_kmeans_patients_by_div_gene/f1_kmeans_PAML_patients_by_PAML_div_genes'.format(project_dir)
clinical_dir = '{}/f0_fran_data_new/data_processed/'.format(project_dir)

# for cohort in ['SG','PAML']:
for cohort in ['PAML']:
    
    # ==== clinical info
    info = pd.read_csv(clinical_dir+os.sep+'{}_clinical.csv'.format(cohort),index_col=0)

    # ==== info of deg clustering
    basename = '{}_pthre5_rk3_ck2'.format(cohort)
    deg_csv = deg_clustering_dir+'/csv/{}.csv'.format(basename)
    deg_list1 = [i.strip() for i in open(deg_clustering_dir+'/txt/{}_genes_cluster1.txt'.format(basename)).readlines()]
    deg_list2 = [i.strip() for i in open(deg_clustering_dir+'/txt/{}_genes_cluster2.txt'.format(basename)).readlines()]
    deg_list3 = [i.strip() for i in open(deg_clustering_dir+'/txt/{}_genes_cluster3.txt'.format(basename)).readlines()]

    patient_list1 = [i.strip() for i in open(patient_clustering_dir+'/txt/{}_patients_cluster1.txt'.format(basename)).readlines()]
    patient_list2 = [i.strip() for i in open(patient_clustering_dir+'/txt/{}_patients_cluster2.txt'.format(basename)).readlines()]

    # ==== re-order the genes followed by down-/div-/up-genes
    reindex = deg_list1+deg_list3+deg_list2

    # position for cluster box
    counts_row = [len(deg_list1),len(deg_list3),len(deg_list2)]

       
    deg_df = pd.read_csv(deg_csv,index_col=0)
    deg_df = deg_df.reindex(reindex)
    # ==== convert the ID into AML_xxx format
    deg_df.columns = info.loc[deg_df.columns]['subject_ID']


    col_colors_dict = {}
    for patient in patient_list1:
        col_colors_dict[info.loc[patient]['subject_ID']]='tab:blue'
    for patient in patient_list2:
        col_colors_dict[info.loc[patient]['subject_ID']]='tab:orange'

    col_colors = []
    for patient in deg_df.columns:
        col_colors.append(col_colors_dict[patient])
        
    for metric in ['euclidean']:
        for method in ['ward',]:
            expr_with_info_plot(deg_df,counts_row,col_colors,mutation_df,cytogenetic_df,outdir,basename,metric,method)

    # for metric in ['euclidean','cosine','correlation']:
    #     for method in ['average','single','weighted','ward','complete','centroid']:
    #         try:
    #             expr_with_info_plot(deg_df,counts_row,col_colors,mutation_df,cytogenetic_df,outdir,basename,metric,method)
    #             print(metric,method,'pass')
    #         except:
    #             print(metric,method,'NOT pass')
    #             pass
    

