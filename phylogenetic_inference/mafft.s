#!/bin/bash
#
# mafft_its2.s
#
# Sergio Alías-Segura
#
#SBATCH -J mafft_its2
#SBATCH -p short
#SBATCH --mem 10G
#SBATCH -t 0-00:10:00 # (D-HH:MM:SS)
#SBATCH -o mafft_its2.%j.out
#SBATCH -e mafft_its2.%j.err

hostname; pwd; date

time mafft --maxiterate 1000 --localpair /mnt/lustre/home/salias/data/secuencias_arboles/Fusarium_ASV_NCBI_outgroup_Cepha.fasta > fus_its2_msa_mafft_linsi.fasta # mafft-linsi
