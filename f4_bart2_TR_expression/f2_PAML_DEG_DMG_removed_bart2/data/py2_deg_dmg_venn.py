import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
#from GenomeData import *
# import AML_modules
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
matplotlib.rcParams["font.family"] = "sans-serif"
import seaborn as sns
sns.set(font_scale=1.5)
sns.set_style("whitegrid", {'axes.grid' : False})
from matplotlib_venn import venn3,venn2
from scipy import stats

#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')


def venn2_plot(a,b,la,lb,figname,cohort,color1='purple',color2='skyblue'):
    # fisher exact test
    total=23251
    s,p = stats.fisher_exact([[len(a.intersection(b)),len(a)],[len(b),total-len(b)]])
    print(cohort,len(a),len(b))
    print('odds ratio={:.2f}, pvalue={:.2e}'.format(s,p))
    
    # plot venn diagram
    plt.figure(figsize=(4,4))
    out = venn2([a,b],set_labels=(la,lb),set_colors=(color1,color2), alpha=0.5)
    
    for text in out.set_labels:
        text.set_fontsize(16)
    for text in out.subset_labels:
        try:
            text.set_fontsize(16)
        except:
            pass
#     plt.title('{}-genes'.format(deg_type))
#     if p<0.05:
#         plt.text(x=.2,y=-.15,s='$p$={:.1e}'.format(p),transform=plt.axes().transAxes)
    plt.savefig(figname,bbox_inches = 'tight',pad_inches=0.1,transparent=True)
    plt.close()
    
    with open(figname+'.txt','w') as outf:
        outf.write('\n'.join(a.intersection(b))+'\n')




# main
    
indir='gene_lists'
outdir = 'f2_deg_dmg_venn'
os.makedirs(outdir,exist_ok=True)

for cohort in ['PAML','SG']:
    deg_list_t = '{}/{}_logFC1_p001_patientGT5_deg_genelist.txt'.format(indir,cohort)
    deg_list_c = '{}/{}_DNAme_genes.txt'.format(indir,cohort)
    
    deg_t = [i.strip() for i in open(deg_list_t).readlines()]
    deg_c = [i.strip() for i in open(deg_list_c).readlines()]
    
    figname = outdir+os.sep+'{}_venn.pdf'.format(cohort)
    venn2_plot(set(deg_t),set(deg_c),'DEG','DMG',figname,cohort)
    




