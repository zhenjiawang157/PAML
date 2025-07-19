import sys,argparse,re
import os,glob
import numpy as np
import pandas as pd
import operator
from collections import Counter
#from GenomeData import *
#import site_align_heatmap
from sklearn.cluster import KMeans
from sklearn.metrics.cluster import adjusted_rand_score

import matplotlib
from matplotlib import gridspec
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.2)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")
sns.set_style("whitegrid", {'axes.grid' : False,'grid.color': 'grey'})
sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'w'})
#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')
# import manu_rerank_genes
# from manu_rerank_genes import label_by_group
# from matplotlib_venn import venn3,venn2
# import AML_modules

def label_by_group(data,groupnum):
    #return new label after k-means clustering
    # newlabels: new labels in a list
    # sublabels: new labels in each group
    kmean = KMeans(groupnum,n_init=20).fit(data)
    label = kmean.labels_  
    labelmap = pd.concat([pd.DataFrame(label,columns=['label']),pd.DataFrame(list(data.index),columns = ['index'])],axis=1)   
    newlabels = []
    sublabels = []
    for i in range(groupnum):
        sub=list(labelmap.loc[labelmap['label']==i]['index'])
        newlabels.extend(sub)
        sublabels.append(sub)
    #print(label,newlabels);exit(0)
    return newlabels,sublabels,label

def return_index_of_consistent_kmeans(df,k,r_times):
    # return the most consistent k-means result from r_times' running
    # label_all is used for the adjusted_rand_score
    index_all,index_sep_all,label_all,ari = [],[],[],{}
    for times in np.arange(r_times): 
        newindex,newindex_sep,newlabel = label_by_group(df,k);#print(newindex,rowlabel);exit(0)
        index_all.append(newindex);#print(newindex_sep)
        index_sep_all.append(newindex_sep)
        label_all.append(newlabel)
    for label_a_index in np.arange(len(label_all)):
        ari[label_a_index] = []
        for label_b_index in np.arange(len(label_all)):
            ari[label_a_index].append(adjusted_rand_score(label_all[label_a_index],label_all[label_b_index]))
            
    for i in sorted(ari.keys(),key=lambda x: np.mean(ari[x]),reverse=True):
        print(df.shape,k,'average ARI: ',np.mean(ari[i]),)
        return index_all[i],index_sep_all[i]


def rerank_by_color(genes,df):

    df = df.loc[genes];#print(df);exit(0)
    df = np.sign(df.abs()).sum(axis=1)
    df = df.sort_values(ascending=False)
    #print(df.index,list(df.index))
    return list(df.index)



def plot_heatmap(plot_df,counts_row,counts_col,basename,outdir):
        
    plt.figure(figsize = (2,3))
    g = sns.heatmap(plot_df,cmap=plt.cm.PiYG_r,cbar=False,vmax=2,vmin=-2,xticklabels=False,yticklabels=False)#,cbar_kws={"orientation": "horizontal"})#"use_gridspec":False,"location":"top"})    
    # plt.setp(g.ax_heatmap.get_xticklabels(),size = 18)# for sns.clustermap
    #plt.setp(g.get_xticklabels(),size = 8)
    plt.ylabel('{} genes'.format(plot_df.shape[0]))
    plt.xlabel('{} patients'.format(plot_df.shape[1]))
    
    for i in np.arange(len(counts_row)+1):
        print(sum(counts_row[:i]))
        plt.axhline(y=sum(counts_row[:i]),xmin=-.00,xmax=1.00,color='k',lw=1,clip_on=False)
    
    for i in np.arange(len(counts_col)+1):
        print(sum(counts_col[:i]))
        plt.axvline(x=sum(counts_col[:i]),ymin=-.00,ymax=1.00,color='k',lw=1,clip_on=False)
    #plt.title(basename)
    plt.savefig(outdir+os.sep+basename+'.png',bbox_inches = 'tight',pad_inches=0.1,transparent=True,dpi=600)
#     plt.show()
    plt.close()

    
def kmeans_clustering_index(df,target_genes,rk,ck,t_times,basename,outdir):

    df_newindex,df_newindex_sep = return_index_of_consistent_kmeans(df,rk,r_times=t_times)
    patient_cluster_genes = df.index.intersection(target_genes)
    df_newcol,df_newcol_sep = return_index_of_consistent_kmeans(np.transpose(df.loc[patient_cluster_genes]),ck,r_times=t_times)    

    # rerank all genes
    df_newindex2,df_newindex_sep2=[],[]
    counts_row,counts_col = [],[]
    for i in np.arange(len(df_newindex_sep)):
        genes = rerank_by_color(df_newindex_sep[i],df)
        df_newindex2=df_newindex2+genes
        df_newindex_sep2.append(genes)
        counts_row.append(len(genes))
    # rerank all patients
    df_newcol2,df_newcol_sep2=[],[]
    for i in np.arange(len(df_newcol_sep)):
        patients = rerank_by_color(df_newcol_sep[i],np.transpose(df))
        df_newcol2=df_newcol2+patients
        df_newcol_sep2.append(patients)
        counts_col.append(len(patients))

    
    #write out the genes lists
    os.makedirs(outdir+os.sep+'txt',exist_ok=True)
    for i in np.arange(len(df_newindex_sep2)):
        with open(outdir+os.sep+'txt'+os.sep+'{}_genes_cluster{}.txt'.format(basename,i+1),'w') as outf:
            outf.write('\n'.join(df_newindex_sep2[i])+'\n')
    for i in np.arange(len(df_newcol_sep2)):
        with open(outdir+os.sep+'txt'+os.sep+'{}_patients_cluster{}.txt'.format(basename,i+1),'w') as outf:
            outf.write('\n'.join(df_newcol_sep2[i])+'\n')
    # genes used for patient cluster
    with open(outdir+os.sep+'txt'+os.sep+'genes_for_patients_cluster.txt','w') as outf:
        outf.write('\n'.join(patient_cluster_genes)+'\n')
            
    return df_newindex2, df_newcol2, counts_row,counts_col

def main():
        
    indir='../f0_deg/f1_cohort_deg/'
    outdir='f1_kmeans_PAML_patients_by_PAML_div_genes'
    os.makedirs(outdir,exist_ok=True)
    os.makedirs(outdir+os.sep+'csv',exist_ok=True)
    
    target_genes = [i.strip() for i in open('../f1_kmeans_patients_by_deg/f1_kmeans_PAML_patients_by_PAML_DEG/txt/PAML_pthre5_rk3_ck2_genes_cluster3.txt').readlines()]
    assert(len(target_genes)<1000)
    
    for patient_thre in[5]:
        for cohort in ['PAML']:
            prename = '{}'.format(cohort)
            logFC = pd.read_csv(indir+os.sep+'{}_logFC1_p001_patientGT{}_deg_logFC.csv'.format(cohort,patient_thre),index_col=0)
            logFC_ori = pd.read_csv(indir+os.sep+'{}_logFC1_p001_patientGT{}_deg_logFC_ori.csv'.format(cohort,patient_thre),index_col=0)
            qvalue = pd.read_csv(indir+os.sep+'{}_logFC1_p001_patientGT{}_deg_qvalue.csv'.format(cohort,patient_thre),index_col=0)

            df = logFC_ori
            t_times = 20
            for rk in [3]:
                for ck in [2]:
                    basename = '{}_pthre{}_rk{}_ck{}'.format(cohort,patient_thre,rk,ck);print('\n\n',basename,'\n')
                    newindex,newcol,counts_row,counts_col = kmeans_clustering_index(df,target_genes,rk,ck,t_times,basename,outdir)
                    plot_df=np.transpose(np.transpose(df).reindex(index=newcol)).reindex(index=newindex)
                    plot_df.to_csv(outdir+os.sep+'csv'+os.sep+basename+'.csv')
                    plot_heatmap(plot_df,counts_row,counts_col,basename,outdir)
            
            
            
    
                    
#                     exit()

        

#     fig,ax = plt.subplots(figsize=(.3,.9))
#     norm = matplotlib.colors.Normalize(vmin=-2, vmax=2)
# #     pal = sns.light_palette('red',as_cmap=True)
#     cb = matplotlib.colorbar.ColorbarBase(ax,cmap=plt.cm.PiYG_r,norm=norm,orientation='vertical')
# #     cb.set_label('log2 FC',rotation=90,labelpad=-32)
#     cb.set_ticks([-2,2])
#     #cb.set_clim(vmax*.15,vmax)
#     ax.tick_params(axis='x',direction='out', length=0, width=1, colors='black')    
#     plt.savefig(outdir+os.sep+'colorbar.png',bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
#     plt.close()





if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', action = 'store', type = str,dest = 'infile', help = 'input file of', metavar = '<file>')
    parser.add_argument('-o','--outfile', action = 'store', type = str,dest = 'outfile', help = 'outfile of', metavar = '<file>')
    #parser.add_argument('-i', '--indir', action = 'store', type = str,dest = 'indir', help = 'input dir of ', metavar = '<dir>')
    #parser.add_argument('-o','--outdir', action = 'store', type = str,dest = 'outdir', help = 'outdir of ,default: current dir', metavar = '<dir>',default='./')
    #parser.add_argument('-s','--species', action = 'store', type = str,dest = 'species', help = 'species used to choose correct chromosome, e.g., hg38 or mm10', metavar = '<str>',required=True)
    

    args = parser.parse_args()
    if(len(sys.argv))<0:
        parser.print_help()
        sys.exit(1)
  
    main()
