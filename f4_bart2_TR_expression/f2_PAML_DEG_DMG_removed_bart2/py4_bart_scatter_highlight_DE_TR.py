import sys,argparse
import os,glob
import numpy as np
import pandas as pd
#from GenomeData import *

import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=15
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=2)
#sns.set_style("whitegrid", {'axes.grid' : False})

#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')


def scatter_plot(df,marked_index,figname,mark_col):
    
    pylist = df.index
    plt.figure(figsize=(3,3))
    for index in df.index:
        xp = list(df.index).index(index)
        values = -1*np.log10(df.loc[index,mark_col])#;print(values)
        plt.scatter(xp,values,color='k',s=6)
    
    max_p = -1*np.log10(df.iloc[0,-1])
    values_reset = max_p*1.05
    # top 3 and TFs in shared top 10
    sorted_marker_ids = np.append(df.index[:10],df.index[:20].intersection(marked_index)[:5]);print(sorted_marker_ids)
    sorted_marker_pos = sorted(set([list(df.index).index(marker) for marker in  sorted_marker_ids]))
    x_pos_m=0
    for marker_id in sorted_marker_pos:
        xp = marker_id #;print(xp,df)
        values = -1*np.log10(df.loc[df.index[marker_id],mark_col])#;print(values)
        plt.scatter(xp,values,color='red',s=15)
        # ==== mark the index label, not overlap with each other
        values_reset = min(values+1,values_reset-max_p*0.09)#;print(values,values_reset)
        # ==== mark the index and plot the arrow separately
        if df.index[marker_id] in marked_index:
            plt.text(xp+120+x_pos_m,values_reset,df.index[marker_id],color='red')
        else:
            plt.text(xp+120+x_pos_m,values_reset,df.index[marker_id])
        
        # ==== plt.arrow(x,y,dx,dy)
        plt.arrow(xp+120+x_pos_m,values_reset+max_p*0.03,-90-x_pos_m,values-values_reset-max_p*0.03,\
                  length_includes_head = True,head_width=max_p*0.02,head_length=35,fc='k',ec='k')
        x_pos_m = x_pos_m + 30

    #plt.title('Predicted co-factors of \nT-ALL gained bindings')
    plt.xlabel('TR Rank')
    plt.ylabel('-log$_{{10}}$ $p$-value')
    plt.axes().set_xticks([0,len(df.index)])
    plt.axes().set_xticklabels([1,len(df.index)],rotation=0, ha='center',color='k')
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,transparent=True,dpi=600)
    plt.show()
    plt.close()




    
# === main
    
    
outdir='f4_bart2_figs'    
os.makedirs(outdir,exist_ok=True)

cohort='PAML'
for group in ['div','up','down']:
    bart_results_file=glob.glob('f1_DMG_removed_DEG_bart2_results/PAML_pthre5_rk3_ck2_genes_*{}_bart_results.txt'.format(group))[0]
    df = pd.read_csv(bart_results_file,sep='\t',index_col=0)
    
    # ttest compr expression between two groups of patients
    expr_ttest1 = pd.read_csv('f3a_PAML_expr_of_topTR_sep_patient_ttest_rel/ttest_PAML_{}_patientG1.csv'.format(group),index_col=0)
    expr_ttest2 = pd.read_csv('f3a_PAML_expr_of_topTR_sep_patient_ttest_rel/ttest_PAML_{}_patientG2.csv'.format(group),index_col=0)
    
    p_thre=0.05
    g1_sig_tfs = expr_ttest1[expr_ttest1['pvalue']<p_thre].index
    g2_sig_tfs = expr_ttest2[expr_ttest2['pvalue']<p_thre].index
    marked_index = g1_sig_tfs.union(g2_sig_tfs).union(['KMT2A'])
    
    # bart scatter plot
    basename = os.path.basename(bart_results_file).split('.txt')[0]
    mark_col = df.columns[-1]
    df = df.sort_values(by=[mark_col],ascending=True)
    figname = outdir+os.sep+basename+'.pdf'
    scatter_plot(df,marked_index,figname,mark_col)


