import sys,argparse
import os,glob
import numpy as np
import pandas as pd
import re,bisect

#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt
#matplotlib.rcParams['font.size']=16
#matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "Bitstream Vera Sans"]
#matplotlib.rcParams["font.family"] = "sans-serif"
#import seaborn as sns
#sns.set(font_scale=2)
#sns.set_style("whitegrid", {'axes.grid' : False})
#sns.set_style("ticks",{'ytick.color': 'k','axes.edgecolor': 'k'})
#sns.despine(offset=0, trim=True)
# import AML_modules


  

def return_deg_info(cohort,logfc_thre,qvalue_thre):
    
    '''This is used to return the deseq2 results and TPM for patient samples '''
    file_dir = '/nv/vol190/zanglab/zw5j/wwork2018/AML_fran/f0_fran_data_new/data_processed'
    
    compr_logFC = pd.read_csv(file_dir+os.sep+'{}_deseq2_log2FC.txt'.format(cohort),sep=',',index_col=0)
    compr_logFC_ori = compr_logFC
    compr_qvalue = pd.read_csv(file_dir+os.sep+'{}_deseq2_qvalue.txt'.format(cohort),sep=',',index_col=0)
        
    # only keep logFC values when qvalue < qvalue_thre and abs(logfc)>logfc_thre
    compr_logFC = compr_logFC[compr_qvalue < qvalue_thre]
    compr_logFC = compr_logFC[(compr_logFC > logfc_thre)|(compr_logFC < -1*logfc_thre)]
    compr_logFC = compr_logFC.fillna(0)
    
    return compr_logFC,compr_logFC_ori,compr_qvalue


def main():
        
    outdir='f1_cohort_deg'
    os.makedirs(outdir,exist_ok=True)
    
    project_dir="/nv/vol190/zanglab/zw5j/wwork2018/AML_fran/"
    clinical_dir = '{}/f0_fran_data_new/data_processed/'.format(project_dir)
    logfc_thre,qvalue_thre = 1,0.01
    for patient_thre in[5]:
        writer = pd.ExcelWriter(outdir+os.sep+'DESeq2_DEG.xlsx')

        for cohort in ['SG','PAML']:
            logFC,logFC_ori,qvalue = return_deg_info(cohort,logfc_thre,qvalue_thre)
            info = pd.read_csv(clinical_dir+os.sep+'{}_clinical.csv'.format(cohort),index_col=0)
            
            sum_logFC = np.sign(logFC.abs()).sum(axis=1);#print(sum_paml)
            deg_index = sum_logFC[sum_logFC > patient_thre].index
            
            logFC.loc[deg_index].to_csv(outdir+os.sep+'{}_logFC1_p001_patientGT{}_deg_logFC.csv'.format(cohort,patient_thre))
            logFC_ori.loc[deg_index].to_csv(outdir+os.sep+'{}_logFC1_p001_patientGT{}_deg_logFC_ori.csv'.format(cohort,patient_thre))
            qvalue.loc[deg_index].to_csv(outdir+os.sep+'{}_logFC1_p001_patientGT{}_deg_qvalue.csv'.format(cohort,patient_thre))
            
            logFC.loc[deg_index].rename(columns=info['subject_ID']).round(5).to_csv(outdir+os.sep+'{}_logFC1_p001_patientGT{}_deg_logFC_subjectID.csv'.format(cohort,patient_thre))
            logFC_ori.loc[deg_index].rename(columns=info['subject_ID']).round(5).to_csv(outdir+os.sep+'{}_logFC1_p001_patientGT{}_deg_logFC_ori_subjectID.csv'.format(cohort,patient_thre))
            qvalue.loc[deg_index].rename(columns=info['subject_ID']).to_csv(outdir+os.sep+'{}_logFC1_p001_patientGT{}_deg_qvalue_subjectID.csv'.format(cohort,patient_thre))
            
            # save to the excel file
            logFC_ori.loc[deg_index].rename(columns=info['subject_ID']).round(5).fillna('NA').to_excel(writer,'{} DEG log2FC'.format(cohort))
            qvalue.loc[deg_index].rename(columns=info['subject_ID']).fillna('NA').to_excel(writer,'{} DEG adj.p'.format(cohort))
              
            with open(outdir+os.sep+'{}_logFC1_p001_patientGT{}_deg_genelist.txt'.format(cohort,patient_thre),'w') as outf:
                outf.write('\n'.join(deg_index)+'\n')
                
        writer.save()
        


        





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
