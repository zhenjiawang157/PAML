import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
#from GenomeData import *
import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.rcParams['font.size']=16
matplotlib.rcParams["font.sans-serif"] = ["Arial"]
import seaborn as sns
sns.set(font_scale=1.1)
sns.set_style("whitegrid", {'axes.grid' : False})
sns.set_style("ticks")
from scipy import stats
#import re,bisect
#plus = re.compile('\+')
#minus = re.compile('\-')


def return_deg(cohort):
    if cohort =="PAML":
        dn_file = 'data/PAML_pthre5_rk3_ck2_gene1_down.txt'
        up_file = 'data/PAML_pthre5_rk3_ck2_gene2_up.txt'  
        div_file = 'data/PAML_pthre5_rk3_ck2_gene3_div.txt' 
    elif cohort =="SG":
        dn_file = 'data/SG_pthre5_rk3_ck2_gene1_down.txt'
        up_file = 'data/SG_pthre5_rk3_ck2_gene2_up.txt'   
        div_file = 'data/SG_pthre5_rk3_ck2_gene3_div.txt'
    divgenes = [i.strip() for i in open(div_file).readlines()]
    upgenes = [i.strip() for i in open(up_file).readlines()]
    dngenes = [i.strip() for i in open(dn_file).readlines()]
    return divgenes,upgenes,dngenes
    
    

def mark_pvalue(compr_pos,positions,box_vals):
    s,p = stats.ttest_ind(box_vals[compr_pos[0]],box_vals[compr_pos[1]],nan_policy='omit')
    y, h, col = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),99.7)*1.01 ,1.05, 'k'
    y2 = np.percentile(np.append(box_vals[compr_pos[0]],box_vals[compr_pos[1]]),0)*0.99
    x1,x2 = positions[compr_pos[0]],positions[compr_pos[1]]
    p_label='*'
    if p<0.05:
        if p<0.001:
            p_label = '**'
        if compr_pos[2] == 't':
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y, y*h, y*h, y], lw=1, c=col)
            plt.text((x1+x2)*.5, y*1.13, '{:.1e}'.format(p), ha='center', va='center', color=col,fontsize=13)
        else:
            plt.plot([x1*1.03, x1*1.03, x2*0.97, x2*0.97], [y2, y2-2, y2-2, y2], lw=1, c=col)
            plt.text((x1+x2)*.5, y2-2, '{:.1e}'.format(p), ha='center', va='center', color=col,fontsize=13)


def rp_compr_boxplot(box_vals,figname,xticklabels,colors,gene_group,factor,term):
    positions = [1,2]
    plt.figure(figsize=(3,3))
    g = plt.boxplot(box_vals,positions=positions,widths = .6,patch_artist=True,\
                boxprops=dict(color='k',facecolor='w',fill=None,lw=1),\
                medianprops=dict(color='grey'),showfliers=False)    
    for patch, color in zip(g['boxes'], colors):
        patch.set_facecolor(color)

    scatter_X = []
    for position_id in np.arange(len(positions)):
        scatter_x = np.random.normal(positions[position_id],0.07,len(box_vals[position_id]))
        plt.scatter(scatter_x,box_vals[position_id],color=colors[position_id],s=10,zorder=0,alpha=0.99,rasterized=True)
    
    for compr_pos in [[0,1,'t']]:
        mark_pvalue(compr_pos,positions,box_vals)

    plt.axes().set_xticklabels(xticklabels,fontsize=16,rotation=30,ha='right')
    plt.ylabel('sqrt(RP)',fontsize=18)
    plt.title('{}\n{}'.format(factor,term),fontsize=16)
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.1,dpi=600,transparent=True)
    plt.show()
    plt.close()





# ==== main section

outdir = 'f2_RP_compr_enrichr_genes_figs'
os.makedirs(outdir,exist_ok=True)

project_dir='/Volumes/zanglab/zw5j/wwork2018/AML_fran'
groups_matches = {'PAML_pthre5_rk3_ck2_gene1_down':'dn',               
                  'PAML_pthre5_rk3_ck2_gene3_div':'div',
                  'PAML_pthre5_rk3_ck2_gene2_up':'up',}

cohort='PAML'
divgenes,upgenes,dngenes = return_deg(cohort)

ttest_df_s = pd.DataFrame()
ttest_df_p = pd.DataFrame()

for group_key in list(groups_matches.keys())[:]:
    gene_group = groups_matches[group_key]
    # gene_group = group_key.split('_')[-1]
    report_file = 'enrichr/{}/KEGG_2019_Human.human.enrichr.reports.txt'.format(group_key)
    report_df = pd.read_csv(report_file,sep='\t',index_col=1)
    for term in report_df.index[:2]:
        enriched_genes = report_df.loc[term,'Genes'].split(';')
    
        rp_files = glob.glob('{}/manuscript_figs/f6_DEG_bart/f4_RP_on_genes/{}_{}/*.txt'.format(project_dir,cohort,gene_group))
        for rp_file in rp_files[:]:
            rp_file_basename = os.path.basename(rp_file).split('.txt')[0];print(rp_file_basename)
            factor = '_'.join(rp_file_basename.split('_')[:2])
            if factor.startswith(('SPI1_GSM2804465','GATA1_GSM722394','TAL1_GSM1038699','RBBP5_GSM831013')):
                rp_df = pd.read_csv(rp_file,sep='\t',header=None,index_col=0)
                rp_df.columns=['rp']
                rp_df = np.sqrt(rp_df)
                div_rps = rp_df.loc[divgenes,'rp'].dropna().values
                dn_rps = rp_df.loc[dngenes,'rp'].dropna().values
                up_rps = rp_df.loc[upgenes,'rp'].dropna().values
                enriched_rps = rp_df.loc[enriched_genes,'rp'].dropna().values
                
                if gene_group=='div': 
                    assert(len(set(enriched_genes).intersection(divgenes))!=0)
                    box_vals = [rp_df['rp'].values,enriched_rps]
                    xticklabels = ['all genes (#{})'.format(rp_df.shape[0]),
                                   'GroupB genes enriched in \n {} (#{})'.format(term, len(enriched_rps))]
                elif gene_group=='up':           
                    assert(len(set(enriched_genes).intersection(upgenes))!=0)
                    box_vals = [rp_df['rp'].values,enriched_rps]
                    xticklabels = ['all genes (#{})'.format(rp_df.shape[0]),
                                   'GroupC genes enriched in \n {} (#{})'.format(term, len(enriched_rps))]
                elif gene_group=='dn':           
                    assert(len(set(enriched_genes).intersection(dngenes))!=0)
                    box_vals = [rp_df['rp'].values,enriched_rps]
                    xticklabels = ['all genes (#{})'.format(rp_df.shape[0]),
                                   'GroupA genes enriched in \n {} (#{})'.format(term, len(enriched_rps))]
                colors = ['k']*2
                figname = '{}/{}_{}_{}.pdf'.format(outdir,gene_group,factor,term)
                rp_compr_boxplot(box_vals,figname,xticklabels,colors,gene_group,factor,term)
                s,p = stats.ttest_ind(box_vals[1],box_vals[0],nan_policy='omit')
                ttest_df_s.loc[factor.split('_')[0],term] = s 
                ttest_df_p.loc[factor.split('_')[0],term] = p
ttest_df_s.to_csv('{}/stats_s.csv'.format(outdir))
ttest_df_p.to_csv('{}/stats_p.csv'.format(outdir))
    
            
            
            
# plot the heatmap
fig = plt.figure(figsize=(3,3))
g1 = sns.heatmap(np.transpose(-1*np.log10(ttest_df_p)),cmap = plt.cm.Reds,xticklabels=True, yticklabels=True,\
                cbar=True,vmin=None,vmax=10,\
                cbar_kws={"ticks":[0,5,10]})#"use_gridspec":False,"location":"top"})    
plt.savefig('{}/stats_p.pdf'.format(outdir),bbox_inches='tight',pad_inches=0.1,dpi=600)
plt.show()
plt.close()

        
        
 




