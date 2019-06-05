#!/bin/bash
#PBS -P te53
#PBS -q express
#PBS -l walltime=04:00:00
#PBS -l mem=16GB
#PBS -l ncpus=1
#PBS -N MCC_impute2
#PBS -m e
#PBS -M u5015730@anu.edu.au
#PBS -j oe
#PBS -o /g/data1a/te53/MitoImpute/logs/

# LOAD THE MODULE
module unload intel-fc intel-cc
module load python/2.7.11
module load intel-fc/16.0.3.210
module load intel-cc/16.0.3.210
module load Rpackages/3.4.3

# SPECIFY REFERENCE PANEL
#REFpanel="ReferencePanel_v5"
REFpanel=${REFpanel}
HAPLOGREP=~/GitCode/MitoImputePrep/haplogrep/2.1.19/haplogrep-2.1.19.jar
echo
echo "REFERENCE PANEL:				${REFpanel}"
echo "SNP CHIP:						${MtPlatforms}"
echo "HAPLOGREP:					${HAPLOGREP}"
echo "MCMC LENGTH:					${mcmc}"
echo "MCMC BURN-IN:					${burn}"
echo "REF HAPLOTYPES (k_hap):		${khap}"
echo "EFFECTIVE POP. SIZE (Ne):		${ne}"
echo

WGS_VCF=/g/data1a/te53/MitoImpute/data/VCF/chrMT_1kg_SNPonly
TYP_VCF=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.vcf.gz

echo IMP_PRE
#IMP_VCF=
#IMP_INFO=
#OUT_PRE=

