#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 48:00:00
#SBATCH -p largemem
#SBATCH -A cphg_cz3d
#SBATCH -o out1_kmeans_DEG_bart2.log


time bart2 geneset -i ../../f1_kmeans_patients_by_deg/f1_kmeans_PAML_patients_by_PAML_DEG/txt/PAML_pthre5_rk3_ck2_genes_cluster1.txt -s hg38 --outdir f1_kmeans_DEG_bart2_results --ofilename PAML_pthre5_rk3_ck2_genes_cluster1_down 
time bart2 geneset -i ../../f1_kmeans_patients_by_deg/f1_kmeans_PAML_patients_by_PAML_DEG/txt/PAML_pthre5_rk3_ck2_genes_cluster2.txt -s hg38 --outdir f1_kmeans_DEG_bart2_results --ofilename PAML_pthre5_rk3_ck2_genes_cluster2_up
time bart2 geneset -i ../../f1_kmeans_patients_by_deg/f1_kmeans_PAML_patients_by_PAML_DEG/txt/PAML_pthre5_rk3_ck2_genes_cluster3.txt -s hg38 --outdir f1_kmeans_DEG_bart2_results --ofilename PAML_pthre5_rk3_ck2_genes_cluster3_div



