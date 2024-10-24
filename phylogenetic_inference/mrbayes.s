#!/bin/bash
#
# mrbayes_its2.s
#
# Sergio Alías-Segura
#
#SBATCH -J mrbayes_its2
#SBATCH -p short
#SBATCH --mem 50G
#SBATCH -t 0-01:00:00 # (D-HH:MM:SS)
#SBATCH -o mrbayes_its2.%j.out
#SBATCH -e mrbayes_its2.%j.err
#SBATCH -n 8  # Request 8 tasks for MrBayes (2 parallel executions, 4 chains per execution)

hostname; pwd; date

time mpiexec -np 8 mb-mpi instrucciones-MB-Fusarium-ASV-NCBI.nex
