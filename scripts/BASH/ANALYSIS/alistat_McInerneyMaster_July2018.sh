#!/bin/bash
#PBS -P gw26
#PBS -q biodev
#PBS -l walltime=168:00:00
#PBS -l mem=12GB
#PBS -l ncpus=1
#PBS -m e
#PBS -M u5015730@anu.edu.au
#PBS -j oe
#PBS -o /g/data1a/te53/MitoImpute/logs/
#PBS -l storage=scratch/te53+gdata/te53
#PBS -N AliStat_McInerneyMaster_July2018


# LOAD THE MODULE FILE
module load AliStat

# SET VARIABLES
FASTA=/g/data1a/te53/MitoImpute/data/FASTA/masters/McInerney_Master_Alignment_July18_2018.fasta
OUT_DIR=/g/data1a/te53/MitoImpute/analyses/AliStat/`basename ${FASTA} .fasta`
OUT_FILE=${OUT_DIR}`basename ${FASTA} .fasta`

if [ ! -d ${OUT_DIR} ]
then
	mkdir -p ${OUT_DIR}
fi

alistat \
	${FASTA} 1 \
	-o ${OUT_FILE} \
	-d


# END !