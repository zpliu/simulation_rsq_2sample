#!/bin/bash
#PBS -l nodes=1:ppn=12
#PBS -l mem=10gb
#PBS -l walltime=00:30:00
#PBS -m abe
#PBS -q medium
#PBS -N hpcf2hpcf2_HEIDIvMIX

cd /home/zpliu/projectIP/simulation_YJ_2sample

thread=12
gwas_thresh=5e-8
heidi_thresh=0.01

### Run flexmix
time Rscript plot_flexmix_unix.R --gwas_thresh $gwas_thresh \
								--thread $thread

### Run HEIDI								
time Rscript HEIDI_outlier_YJ_unix.R --gwas_thresh $gwas_thresh \
									--heidi_thresh $heidi_thresh \
									--thread $thread

### Run comparison
Rscript plot_HEIDIvMIX.R

