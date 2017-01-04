#!/bin/bash
#PBS -S /bin/bash
#PBS -P cegs
#PBS -q new_nodes
#PBS -J 1-42
#PBS -l select=1:ncpus=4:mem=40gb:local_disk=30G
#PBS -l walltime=48:00:00
#PBS -o /home/GMI/rahul.pisupati/logs/logs.bsmeth
#PBS -e /home/GMI/rahul.pisupati/logs/logs.bsmeth


mkdir $LOCAL_DISK/tmp
export TMPDIR="$LOCAL_DISK/tmp/"
#export TMPDIR="/lustre/scratch/users/rahul.pisupati/tempFiles/"

cd $PBS_O_WORKDIR

scriptFol="/home/GMI/rahul.pisupati/MyScripts/pyBsHap/scripts/"

RefSeq="/home/GMI/rahul.pisupati/TAiR10_ARABIDOPSIS/TAIR10_wholeGenome.fasta"

#inFile=`ls *no_clonal.bam | head -n $PBS_ARRAY_INDEX | tail -n 1 | cut -f1 -d "_"`
inFile=`ls *modified.bam | head -n $PBS_ARRAY_INDEX | tail -n 1`
inID=`echo $inFile | cut -f1 -d "_"`

outFol="methReads"
mkdir -p $outFol/logs

module load pysam/0.9.1.4-foss-2016a-Python-2.7.11
module load pyfaidx/0.4.8.1-foss-2016a-Python-2.7.11
module load numpy/1.10.4-foss-2016a-Python-2.7.11

context=(
  'CN'
  'CG'
  'CHG'
  'CHH'
)

for (( i=0; i<4; i=i+1 ));do
  python $scriptFol/potatoskin_getMethWindows.py getmeth -i ${inFile} -r $RefSeq -o $outFol/meths.${inID}.${context[$i]}.json -c ${context[$i]} -v > $outFol/logs/logs.meths.${inID}.${context[$i]} 2>&1 &
done
wait
