import sys,argparse,re
import os,glob
import numpy as np
import pandas as pd
import operator




# == main ==
outdir='data'
os.makedirs(outdir,exist_ok=True)

# read and save beat aml expression data
# project_dir="/nv/vol190/zanglab/zw5j/wwork2018/AML_fran"
project_dir="/Volumes/zanglab/zw5j/wwork2018/AML_fran"
beatAML_dir='{}/manuscript_figs/f11_clinical_TCGA_and_beatAML/beatAML/'.format(project_dir)
# beatAML_RPKM = pd.read_excel('{}/beatAML_RPKM.xlsx'.format(beatAML_dir),index_col=0)
# beatAML_CPM = pd.read_excel('{}/beatAML_CPM.xlsx'.format(beatAML_dir),index_col=0)
# beatAML_RPKM.round(2).to_csv('{}/beatAML_RPKM.csv'.format(beatAML_dir))
# beatAML_CPM.round(2).to_csv('{}/beatAML_CPM.csv'.format(beatAML_dir))
beatAML_RPKM = pd.read_csv('{}/beatAML_RPKM.csv'.format(beatAML_dir),index_col=0)
beatAML_CPM = pd.read_csv('{}/beatAML_CPM.csv'.format(beatAML_dir),index_col=0)

beatAML_clinical_df = pd.read_excel('{}/beatAML_suppl_ClinicalInformation_7+3Tx.xlsx'.format(beatAML_dir),index_col=0,sheet_name='7+3') 
beatAML_clinical_df = beatAML_clinical_df[['overallSurvival','vitalStatus']].dropna()
# ==== check the patients with both clinical data and expression data
shared_patients = beatAML_clinical_df.index.intersection(beatAML_RPKM.columns)
beatAML_RPKM = beatAML_RPKM[shared_patients]
beatAML_CPM = beatAML_CPM[shared_patients]
beatAML_clinical_df = beatAML_clinical_df.loc[shared_patients]
beatAML_clinical_df.loc[beatAML_clinical_df['vitalStatus']=='Dead','vitalStatus']=1
beatAML_clinical_df.loc[beatAML_clinical_df['vitalStatus']=='Alive','vitalStatus']=0
beatAML_clinical_df.columns = ['time','status']
beatAML_clinical_df.to_csv(outdir+os.sep+'clinical_data.csv')

# # read bart results
bart_results = '{}/updated_202010/f4_bart2_TR_expression/f2_PAML_DEG_DMG_removed_bart2/f1_DMG_removed_DEG_bart2_results/PAML_pthre5_rk3_ck2_genes_cluster3_div_bart_results.txt'.format(project_dir)
tfs_df = pd.read_csv(bart_results,sep='\t',index_col=0)
p_thre = 3
top_tfs = tfs_df[tfs_df['irwin_hall_pvalue']<np.power(.1,p_thre)].index;print(len(top_tfs))
tf_RPKM = beatAML_RPKM.loc[top_tfs].dropna()
tf_RPKM = np.transpose(tf_RPKM)
tf_RPKM.to_csv(outdir+os.sep+'top31TF_RPKM.csv')
tf_CPM = beatAML_CPM.loc[top_tfs].dropna()
tf_CPM = np.transpose(tf_CPM)
tf_CPM.to_csv(outdir+os.sep+'top31TF_CPM.csv')
