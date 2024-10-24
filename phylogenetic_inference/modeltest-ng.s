#!/bin/bash
#
# modeltest-ng.s
#
# Sergio Alías-Segura
#
#SBATCH -J modeltest-ng
#SBATCH -p short
#SBATCH --mem 5G
#SBATCH -t 0-00:10:00 # (D-HH:MM:SS)
#SBATCH -o modeltest-ng.%j.out
#SBATCH -e modeltest-ng.%j.err

hostname; pwd; date

time modeltest-ng --datatype nt \
  --input edited_fus_its2_msa_mafft_linsi.fasta \
  --output modeltest_results \
  --rngseed 1234 \
  -t ml \
  --model-het uigf \
  --template mrbayes \
