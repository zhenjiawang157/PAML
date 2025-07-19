import sys,argparse
import os,glob,re
import numpy as np
import pandas as pd
#from GenomeData import *
import matplotlib



# mutation changes of cohort II
data_prename = 'mutation_DvsR_cohortII'
mutation_df = pd.read_excel('41375_2021_1153_MOESM4_ESM.xlsx')
mutation_df['mutation_score'] = mutation_df.Vaf_Rela - mutation_df.Vaf_Diag
mutation_df.to_csv('{}.csv'.format(data_prename),index=False)

out_df = pd.DataFrame()
for ii in mutation_df.index:
    patient = mutation_df.loc[ii,'SAMPLE_NAME'] 
    gene = mutation_df.loc[ii,'GENE'] #;print(ii,gene)
    mutation_score = mutation_df.loc[ii,'mutation_score']
    if mutation_score > 0.05:
        out_df.loc[patient,gene] = 1
    elif mutation_score < -0.05:
        out_df.loc[patient,gene] = -1
    else:
        out_df.loc[patient,gene] = 0
out_df.to_csv('{}_matrix.csv'.format(data_prename))
            

# mutation changes of cohort I
data_prename = 'mutation_DvsR_cohortI'
mutation_df = pd.read_excel('41375_2021_1153_MOESM7_ESM.xlsx')
mutation_df['mutation_score'] = mutation_df['VAF at relapse'] - mutation_df['VAF at diagnosis']
mutation_df.to_csv('{}.csv'.format(data_prename),index=False)

out_df = pd.DataFrame()
for ii in mutation_df.index:
    patient = mutation_df.loc[ii,'Study Subject ID'] 
    gene = mutation_df.loc[ii,'Gene'] #;print(ii,gene)
    mutation_score = mutation_df.loc[ii,'mutation_score']
    if mutation_score > 0.05:
        out_df.loc[patient,gene] = 1
    elif mutation_score < -0.05:
        out_df.loc[patient,gene] = -1
    else:
        out_df.loc[patient,gene] = 0
out_df.to_csv('{}_matrix.csv'.format(data_prename))
            


# CNV changes of cohort I
data_prename = 'CNV_DvsR_cohortI'
df = pd.read_excel('41375_2021_1153_MOESM6_ESM.xlsx')

out_df = pd.DataFrame()
for ii in df.index:
    patient = df.loc[ii,'Individual'] 
    event = df.loc[ii,'Event'] 
    chrom = df.loc[ii,'Chr'] 
    start = df.loc[ii,'Start'] 
    end = df.loc[ii,'End'] 
    col = '{}_{}_{}_{}'.format(event,chrom,start,end)
    stage = df.loc[ii,'Disease stage']
    if stage =='relapse':
        out_df.loc[patient,col] = 1
    elif stage =='diagnosis':
        out_df.loc[patient,col] = -1
    elif stage =='both':
        out_df.loc[patient,col] = 0
out_df.to_csv('{}_matrix.csv'.format(data_prename))



# genetic and epigenetic evolution models 
data_prename = 'genetic_evolution'
df = pd.read_excel('progressionModelsAnnotation2021.xlsx')
df.columns = ['SubjectID','ID','epiAnn','epiallele.shift',
              'epigenetic_cluster','genetic_cluster','cohortI_pc','cohortII_pc']

df.loc[df.genetic_cluster=='Clonal changes','genetic_cluster_new']=1
df.loc[df.genetic_cluster=='Stable','genetic_cluster_new']=2
df.loc[df.genetic_cluster=='Subclonal changes','genetic_cluster_new']=3

df1 = df[df.cohortI_pc.notnull()]
df2 = df[df.cohortII_pc.notnull()]
df1.to_csv('evolution_corhotI.csv'.format(data_prename),index=False)
df2.to_csv('evolution_corhotII.csv'.format(data_prename),index=False)




