import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
#from GenomeData import *
# import AML_modules
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=12
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.1)
sns.set_style("whitegrid", {'axes.grid' : False})
from matplotlib_venn import venn3,venn2
from scipy import stats

#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')


def venn2_plot(a,b,deg_type,figname,outdir,stat_data,color1='purple',color2='skyblue'):
    # fisher exact test
    total=23251
    s,p = stats.fisher_exact([[len(a.intersection(b)),len(a)],[len(b),total-len(b)]])
    stat_data.loc[deg_type,'fisher_exact_s'] = s
    stat_data.loc[deg_type,'fisher_exact_p'] = p
    # print(deg_type,len(a),len(b))
    # print('odds ratio={:.2f}, pvalue={:.2e}'.format(s,p))
    
    # plot venn diagram
    plt.figure(figsize=(2.6,2.6))
    out = venn2([a,b],set_labels=(deg_type,'MLL targets'),set_colors=(color1,color2), alpha=0.5)

    # reset the text
    label= out.get_label_by_id('01').get_text()
    out.get_label_by_id('01').set_text('    {}'.format(label))
    
    for text in out.set_labels:
        text.set_fontsize(12)
    for text in out.subset_labels:
        try:
            text.set_fontsize(12)
        except:
            pass

    # plt.title('{}-genes'.format(deg_type))
    if p<0.05:
        plt.text(x=.2,y=-.15,s='$p$={:.1e}'.format(p),transform=plt.axes().transAxes)
    plt.savefig(figname+'.pdf',bbox_inches = 'tight',pad_inches=0.1,transparent=True)
    plt.close()
    
    with open(figname+'.txt','w') as outf:
        outf.write('\n'.join(a.intersection(b))+'\n')

    return stat_data


# ==== main
outdir = 'deg_venn_MLL_targets'
os.makedirs(outdir,exist_ok=True)

# group ABC genes    
project_dir="/nv/vol190/zanglab/zw5j/wwork2018/AML_fran"
project_dir="/Volumes/zanglab/zw5j/wwork2018/AML_fran"
deg_clustering_dir = '{}/updated_202010/f1_kmeans_patients_by_deg/f1_kmeans_PAML_patients_by_PAML_DEG'.format(project_dir)
basename = 'PAML_pthre5_rk3_ck2'
deg_csv = deg_clustering_dir+'/csv/{}.csv'.format(basename)
deg_list1 = [i.strip() for i in open(deg_clustering_dir+'/txt/{}_genes_cluster1.txt'.format(basename)).readlines()]
deg_list2 = [i.strip() for i in open(deg_clustering_dir+'/txt/{}_genes_cluster2.txt'.format(basename)).readlines()]
deg_list3 = [i.strip() for i in open(deg_clustering_dir+'/txt/{}_genes_cluster3.txt'.format(basename)).readlines()]
deg_dict = {'groupA':deg_list1,
            'groupB':deg_list3,
            'groupC':deg_list2}

# === MLL target genes
mll_targets = pd.read_excel('MLLTc_Targetgenes.xlsx',index_col=0)
stat_data = pd.DataFrame()

for deg_type in ['groupA','groupB','groupC']:
    deg_t = deg_dict[deg_type]
    deg_c = mll_targets.index
    figname = outdir+os.sep+'{}_venn_MLL_targets'.format(deg_type)
    stat_data = venn2_plot(set(deg_t),set(deg_c),deg_type,figname,outdir,stat_data)

stat_data.to_csv(outdir+os.sep+'stat_data.csv')


