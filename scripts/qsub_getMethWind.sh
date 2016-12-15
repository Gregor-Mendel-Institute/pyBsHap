#!/bin/bash
#PBS -S /bin/bash
#PBS -P cegs
#PBS -q new_nodes
#PBS -l select=1:ncpus=1:mem=10gb:local_disk=10G
#PBS -l walltime=48:00:00
#PBS -o /home/GMI/rahul.pisupati/logs/logs.bsmeth
#PBS -e /home/GMI/rahul.pisupati/logs/logs.bsmeth


mkdir $LOCAL_DISK/tmp
export TMPDIR="$LOCAL_DISK/tmp/"
#export TMPDIR="/lustre/scratch/users/rahul.pisupati/tempFiles/"

cd $PBS_O_WORKDIR

scriptFol="/home/GMI/rahul.pisupati/MyScripts/pyBsHap/scripts/"

RefSeq="/home/GMI/rahul.pisupati/TAiR10_ARABIDOPSIS/TAIR10_wholeGenome.fasta"

inFile=`ls *no_clonal.bam | head -n 1 | cut -f1 -d "_"`

python $scriptFol/potatoskin_getMethWindows.py getmeth -i ${inFile}_processed_reads_no_clonal.bam -r $RefSeq -o meths.${inFile}.json -v > logs.meths.${inFile} 2>&1
