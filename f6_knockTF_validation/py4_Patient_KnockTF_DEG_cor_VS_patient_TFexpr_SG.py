import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
import seaborn as sns
sns.set(font_scale=1.5)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
matplotlib.rcParams['mathtext.fontset'] = 'custom'
matplotlib.rcParams["mathtext.rm"] = "Arial"
from scipy import stats



def return_patient_info(cohort,project_dir):
    '''This is used to return the deseq2 results and TPM for patient samples '''
    file_dir = '{}/f0_fran_data_new/data_processed'.format(project_dir)
    compr_logFC = pd.read_csv(file_dir+os.sep+'{}_deseq2_log2FC.txt'.format(cohort),sep=',',index_col=0)
    compr_qvalue = pd.read_csv(file_dir+os.sep+'{}_deseq2_qvalue.txt'.format(cohort),sep=',',index_col=0)
    return compr_logFC,compr_qvalue


def scatter_plot(out_df,outdir,cohort,dataset,tf,patient_down,patient_up):

    plt.figure(figsize=(2.6,2.6))
    x,y = out_df['TFexpr'],out_df['coefficient']
    # t = plt.scatter(x,y,s=22,c='k',rasterized=True)
    # x2,y2 = out_df[out_df['pvalue']<0.01]['TFexpr'],out_df[out_df['pvalue']<0.01]['coefficient']
    for ii in out_df.index:
        xi,yi,pi = out_df.at[ii,'TFexpr'],out_df.at[ii,'coefficient'],out_df.at[ii,'pvalue']
        if ii in patient_down and pi<0.05:
            plt.scatter(xi,yi,s=22,c='tab:blue',rasterized=True)
        elif ii in patient_up and pi<0.05:
            plt.scatter(xi,yi,s=22,c='tab:orange',rasterized=True)
        else :
            plt.scatter(xi,yi,s=22,c='grey',rasterized=True)
    # regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)       
    x_sort = np.sort(x)
    plt.plot(x_sort,x_sort*slope+intercept,c = 'grey',ls='--',lw=.6)
    if dataset=='GSE87055':
        plt.text(.1,.87,'$R^2$ = {:.2f}'.format(r_value**2),fontsize=13,transform=plt.axes().transAxes)
    else:
        plt.text(.56,.87,'$R^2$ = {:.2f}'.format(r_value**2),fontsize=13,transform=plt.axes().transAxes)

    plt.axhline(y=0,color='k',lw=1.2,ls='--')
    plt.axvline(x=0,color='k',lw=1.2,ls='--')
    plt.title('{}'.format(cohort),fontsize=14)
    plt.xlabel('{} differential expression\n in patients'.format(tf),fontsize=17)
    plt.ylabel('Correlation between \npatient and knockTF data',fontsize=17)
    figname = '{}/{}_{}_{}.pdf'.format(outdir,cohort,tf,dataset)
    plt.show()
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.close()
    return r_value,p_value


def scatter_plot_pvalue(out_df,outdir,cohort,dataset,tf):

    plt.figure(figsize=(2.6,2.6))
    x,y = out_df['TFexpr'],-1*np.log10(out_df['pvalue'])
    t = plt.scatter(x,y,s=22,c='k',rasterized=True)
    x2,y2 = out_df[out_df['pvalue']<0.01]['TFexpr'],-1*np.log10(out_df[out_df['pvalue']<0.01]['pvalue'])
    t = plt.scatter(x2,y2,s=22,c='red',rasterized=True)

    plt.axhline(y=-1*np.log10(0.01),color='red',lw=1.2,ls='--')
    plt.axvline(x=0,color='k',lw=1.2,ls='--')
#     plt.title('{} {} \n predicted from {} groups'.format(cohort,tf,group),fontsize=14)
    plt.xlabel('{} differential expression\n in patients'.format(tf),fontsize=17)
    plt.ylabel('-log$_{{10}}$($p$value)',fontsize=17)
    figname = '{}/{}_{}_{}_pvalue.pdf'.format(outdir,cohort,tf,dataset)
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.close()
    return 
   

# ==== main section

outdir='f4_Patient_KnockTF_DEG_cor_VS_patient_TFexpr_SG'
os.makedirs(outdir,exist_ok=True)
project_dir="/nv/vol190/zanglab/zw5j/wwork2018/AML_fran/"
# project_dir="/Volumes/zanglab/zw5j/wwork2018/AML_fran/"
# deg_clustering_dir = '{}/updated_202010/f1_kmeans_patients_by_deg/f1_kmeans_PAML_patients_by_PAML_DEG'.format(project_dir)
# patient_clustering_dir = '{}/updated_202010/f2_kmeans_patients_by_div_gene/f3_kmeans_SG_patients_by_PAML_div_SG_DEG'.format(project_dir)
deg_clustering_dir = '{}/updated_202010/f1_kmeans_patients_by_deg/f2_kmeans_SG_patients_by_PAML_DEG'.format(project_dir)
patient_clustering_dir = '{}/updated_202010/f2_kmeans_patients_by_div_gene/f3_kmeans_SG_patients_by_PAML_div_SG_DEG'.format(project_dir)

knockData = pd.read_excel('{}/manuscript_figs/f8_tf_RP_vs_knockTF_expr/f0_knockTF_data/KnockTF-Browse_selected.xlsx'.format(project_dir),index_col=0)
clinical_dir = '{}/f0_fran_data_new/data_processed/'.format(project_dir)



out_df_sum = pd.DataFrame()               
for dataset in knockData.index[:]:
#for dataset in ['DataSet_01_217']:
    tf = knockData.loc[dataset,'TF']
    if tf not in ['TAL1','MYB','ETS1','GLI2','SPI1','GATA1']:
#     if tf not in ['TAL1']:
        continue
    tissue = knockData.loc[dataset,'TissueType']
    celltype = knockData.loc[dataset,'BiosampleName']
    kd = knockData.loc[dataset,'Knock-Method']
    expr_file='{}/manuscript_figs/f8_tf_RP_vs_knockTF_expr/f0_knockTF_data/differential_datasets/{}.txt'.format(project_dir,dataset)
    if dataset=='GSE74999':
        knockTF_expr = pd.read_excel('{}/manuscript_figs/f8_tf_RP_vs_knockTF_expr/f0_knockTF_data/differential_datasets/GSE74999.xlsx'.format(project_dir),index_col=0)
        knockTF_expr['Log2FC']=knockTF_expr['log2.fold_change.']
    elif dataset=='GSE74999_STAR' or dataset=='GSE74999_K562_PU1_OE':
        knockTF_expr = pd.read_csv('{}/manuscript_figs/f8_tf_RP_vs_knockTF_expr/f0_knockTF_data/differential_datasets/{}.csv'.format(project_dir,dataset),index_col=0)
        knockTF_expr['Log2FC']=knockTF_expr['log2FoldChange'].dropna();print(knockTF_expr)
    elif dataset=='GSE87055':
        knockTF_expr = pd.read_excel('{}/manuscript_figs/f8_tf_RP_vs_knockTF_expr/f0_knockTF_data/differential_datasets/GSE87055.xlsx'.format(project_dir),index_col=1)
        knockTF_expr['Log2FC']=knockTF_expr['logFC']
    else:
        knockTF_expr = pd.read_csv(expr_file,index_col=2,sep='\t')
     
    for cohort in ["SG"]:
        out_df = pd.DataFrame()
        # obtain deg
        basename = '{}_pthre5_rk3_ck2'.format(cohort)
        deg_list1 = [i.strip() for i in open(deg_clustering_dir+'/txt/{}_genes_cluster1.txt'.format(basename)).readlines()]
        deg_list2 = [i.strip() for i in open(deg_clustering_dir+'/txt/{}_genes_cluster2.txt'.format(basename)).readlines()]
        deg_list3 = [i.strip() for i in open(deg_clustering_dir+'/txt/{}_genes_cluster3.txt'.format(basename)).readlines()]
        # ==== re-order the genes followed by down-/div-/up-genes
        deg = deg_list1+deg_list2+deg_list3
        
        patient_down = [i.strip() for i in open(patient_clustering_dir+'/txt/{}_patients_cluster1.txt'.format(basename)).readlines()]
        patient_up = [i.strip() for i in open(patient_clustering_dir+'/txt/{}_patients_cluster2.txt'.format(basename)).readlines()]


        # logFC of patients
        RvsD_logFC,RvsD_qvalue = return_patient_info(cohort,project_dir)
        
        # target genes used for coefficient
        target_genes = RvsD_logFC.index.intersection(knockTF_expr.index)
        target_genes = knockTF_expr.index.intersection(deg)
        
        info = pd.read_csv(clinical_dir+os.sep+'{}_clinical.csv'.format(cohort),index_col=0)
        for patient in RvsD_logFC.columns:
            a = RvsD_logFC[patient][target_genes].fillna(0)
            try:
                b = knockTF_expr['Log2FC'][target_genes].fillna(0) # for knockTF data
            except:
                b = knockTF_expr['Log2FC'].loc[target_genes].drop_duplicates().fillna(0)
            pearsonr_r,pearsonr_p = stats.pearsonr(a,b)  
            out_df.loc[patient,'ID'] = info.loc[patient]['subject_ID']
            out_df.loc[patient,'coefficient'] = pearsonr_r
            out_df.loc[patient,'pvalue'] = pearsonr_p
            out_df.loc[patient,'TFexpr'] = RvsD_logFC.loc[tf][patient]
        out_df.to_csv('{}/_{}_{}_{}.csv'.format(outdir,cohort,tf,dataset))
        r_value,p_value = scatter_plot(out_df,outdir,cohort,dataset,tf,patient_down,patient_up)
#         scatter_plot_pvalue(out_df,outdir,cohort,dataset,tf)
        
        out_df_sum.loc['{}_{}'.format(cohort,dataset),'tf']=tf
        out_df_sum.loc['{}_{}'.format(cohort,dataset),'r']=r_value
        out_df_sum.loc['{}_{}'.format(cohort,dataset),'p']=p_value
#         exit()
        
out_df_sum.to_csv('{}/summary.csv'.format(outdir))
        
        




