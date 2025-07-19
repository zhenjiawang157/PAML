import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
#from GenomeData import *




for cohort in ["PAML"]:
    dmg_file = '/Volumes/zanglab/zw5j/wwork2018/AML_fran/manuscript_figs/f2_DNA_methylation/f1_DNA_methylated_genes/{}_DNAme_genes.txt'.format(cohort)
    dmgs = [i.strip() for i in open(dmg_file).readlines()]
    
    deg_dir='/Volumes/zanglab/zw5j/wwork2018/AML_fran/updated_202010/f1_kmeans_patients_by_deg/f1_kmeans_PAML_patients_by_PAML_DEG'
    deg_lists = glob.glob('{}/txt/{}*genes*.txt'.format(deg_dir,cohort))
    for deg_list in deg_lists:
        basename = os.path.basename(deg_list).split('.txt')[0]
        degs = [i.strip() for i in open(deg_list).readlines()]
        degs_new = [i for i in degs if i not in dmgs]
        print(len(degs),len(degs_new))
        
        with open('PAML_DMG_removed_DEGs/DMG_removed_{}.txt'.format(basename),'w') as outf:
            outf.write('\n'.join(degs_new)+'\n')
