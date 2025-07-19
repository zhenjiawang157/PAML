import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
import scipy
#from GenomeData import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=9
import seaborn as sns
sns.set(font_scale=1)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")
from matplotlib import gridspec
#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
from  scipy.cluster.hierarchy import fcluster



def plot_correlation_figs(cohort,outdir,project_dir,deg_clustering_dir,patient_clustering_dir,flag):
    
    # ==== info of deg clustering
    basename = '{}_pthre5_rk3_ck2'.format(cohort)
    deg_csv = patient_clustering_dir+'/csv/{}.csv'.format(basename)
    deg_df = pd.read_csv(deg_csv,index_col=0)
    # deg_df = np.sqrt(deg_df)
    
    # if use all DEG or div genes to calculate the pairwise correlation
    if flag=='div_genes':
        div_genes = [i.strip() for i in open(deg_clustering_dir+'/txt/{}_genes_cluster3.txt'.format(basename)).readlines()]
        deg_df = deg_df.loc[div_genes]
        
        
    patient_g1 = [i.strip() for i in open(patient_clustering_dir+'/txt/{}_patients_cluster2.txt'.format(basename)).readlines()]
    patient_g2 = [i.strip() for i in open(patient_clustering_dir+'/txt/{}_patients_cluster1.txt'.format(basename)).readlines()]
    
    intra_distances = []
    inter_distances = []
    
    # ==== calculate inter patient dis
    for i in patient_g1:
        for j in patient_g2:
            u = deg_df[i].values
            v = deg_df[j].values
            euclidean_distance = scipy.spatial.distance.euclidean(u, v)
            inter_distances.append(euclidean_distance)
    
    # ==== calculate intra patient dis
    for patient_list in [patient_g1,patient_g2]:
        for i in np.arange(len(patient_list)-1):
            for j in np.arange(i+1,len(patient_list)):
                u = deg_df[patient_list[i]].values
                v = deg_df[patient_list[j]].values
                euclidean_distance = scipy.spatial.distance.euclidean(u, v)
                intra_distances.append(euclidean_distance)
    
    s,p = scipy.stats.ttest_ind(inter_distances,intra_distances)
    print(s,p)
    print(len(inter_distances))
    print(len(intra_distances))
    
    # ==== plot the distribution
    
    fig = plt.figure(figsize=(4,3))
    sns.distplot(intra_distances,label='intra-group',hist=False,color='r')
    sns.distplot(inter_distances,label='inter-group',hist=False,color='grey')
    # plt.axes().set_xticks([df1.columns[0],df1.columns[-1]])
    # plt.axes().set_xticklabels(['-1kb','1kb'],rotation=30)
    plt.ylabel('PDF',fontsize=16)
    plt.xlabel('Euclidean distance',fontsize=16)
    plt.legend(fontsize=14,borderaxespad=0.2,labelspacing=.2,handletextpad=0.2,handlelength=1,loc="upper right",frameon=False)
    plt.text(.65,.5,'p={:.1e}'.format(p),fontsize=14,transform=fig.transFigure)
#     plt.show()
    plt.title(flag)
    plt.savefig(outdir+os.sep+basename+'_{}.pdf'.format(flag),bbox_inches='tight',pad_inches=0.1,dpi=600)
    plt.close()




## ==== main section
    
outdir = 'f3_PAML_patient_correlation'
os.makedirs(outdir,exist_ok=True)

project_dir="/nv/vol190/zanglab/zw5j/wwork2018/AML_fran"
# project_dir="/Volumes/zanglab/zw5j/wwork2018/AML_fran"
cohort='PAML'
# dir of DEG and clinical info
# deg_clustering_dir = '{}/updated_202010/f1_kmeans_patients_by_deg/f1_kmeans_PAML_patients_by_PAML_DEG'.format(project_dir)
patient_clustering_dir = '{}/updated_202010/f2_kmeans_patients_by_div_gene/f1_kmeans_PAML_patients_by_PAML_div_genes'.format(project_dir)        
deg_clustering_dir = '{}/updated_202010/f1_kmeans_patients_by_deg/f1_kmeans_PAML_patients_by_PAML_DEG'.format(project_dir)

plot_correlation_figs(cohort,outdir,project_dir,deg_clustering_dir,patient_clustering_dir,'div_genes')
plot_correlation_figs(cohort,outdir,project_dir,deg_clustering_dir,patient_clustering_dir,'all_degs')


