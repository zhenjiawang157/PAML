#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 48:00:00
#SBATCH -p largemem
#SBATCH -A cphg_cz3d
#SBATCH -o out2_re_kmeans_DEG_bart2.log


time bart2 geneset -i ../../f2_kmeans_patients_by_div_gene/f1_kmeans_PAML_patients_by_PAML_div_genes/txt/PAML_pthre5_rk3_ck2_genes_cluster1.txt -s hg38 --outdir f2_re_kmeans_DEG_bart2_results --ofilename PAML_pthre5_rk3_ck2_genes_cluster1_up 
time bart2 geneset -i ../../f2_kmeans_patients_by_div_gene/f1_kmeans_PAML_patients_by_PAML_div_genes/txt/PAML_pthre5_rk3_ck2_genes_cluster2.txt -s hg38 --outdir f2_re_kmeans_DEG_bart2_results --ofilename PAML_pthre5_rk3_ck2_genes_cluster2_down
time bart2 geneset -i ../../f2_kmeans_patients_by_div_gene/f1_kmeans_PAML_patients_by_PAML_div_genes/txt/PAML_pthre5_rk3_ck2_genes_cluster3.txt -s hg38 --outdir f2_re_kmeans_DEG_bart2_results --ofilename PAML_pthre5_rk3_ck2_genes_cluster3_div



