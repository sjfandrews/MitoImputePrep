#!/bin/bash
#PBS -P te53
#PBS -q express
#PBS -l walltime=00:10:00
#PBS -l mem=32GB
#PBS -l ncpus=1
#PBS -N impute_SNPchip_1kGP
#PBS -m e
#PBS -M u5015730@anu.edu.au
#PBS -j oe
#PBS -o /g/data1a/te53/MitoImpute/logs/

# LOAD THE MODULE
module unload intel-fc intel-cc
module load intel-fc/16.0.3.210
module load intel-cc/16.0.3.210
module load Rpackages/3.4.3
module load bcftools/1.8
module load plink/1.9
module load impute2/2.3.2
module load vt
echo
echo "LOADED R v3.4.3"
echo "LOADED bcftools v1.8"
echo "LOADED plink v1.9"
echo "LOADED IMPUTE2 v2.3.2"

# SPECIFY REFERENCE PANEL
REFpanel="ReferencePanel_v4"
echo
echo "REFERENCE PANEL: ${REFpanel}"

# CHECK FOR OR CREATE THE DECOMPOSED 1,000 GENOMES VCF FILE
norm_vcf=/g/data1a/te53/MitoImpute/data/VCF/chrMT_1kg_norm.vcf.gz
decom_vcf=/g/data1a/te53/MitoImpute/data/VCF/chrMT_1kg_norm_decomposed.vcf.gz
vcf_1kg=/g/data1a/te53/MitoImpute/data/VCF/chrMT_1kg_norm_decomposed_firstAlt.vcf.gz
plink_1kg=/g/data1a/te53/MitoImpute/data/PLINK/chrMT_1kg_norm_decomposed_firstAlt
orig_vcf=/g/data1a/te53/haploco/data/originals/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz
samps_1kg=/g/data1a/te53/MitoImpute/metadata/SampleList1kg.txt
sex_1kg=/g/data1a/te53/MitoImpute/metadata/SampleList1kg_sex.txt
ref_fasta=/g/data1a/te53/MitoImpute/data/FASTA/rCRS.fasta
if [ -f ${vcf_1kg} ]
then
	echo
	echo "${vcf_1kg} EXISTS ... PASSING"
	bcftools query -l ${vcf_1kg} > ${samps_1kg}
	Rscript ~/GitCode/MitoImputePrep/scripts/R/assign_sex_label.R ${samps_1kg} ${sex_1kg}
	cp ${vcf_1kg} ~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/
	cp ${sex_1kg} ~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/
else
	echo
	echo "${vcf_1kg} NOT FOUND ... CREATE THE DECOMPOSED 1,000 GENOMES VCF FILE"
	bcftools norm -f ${ref_fasta} | -m - ${orig_vcf} | bcftools view -V indels,mnps | bcftools norm -m + | bcftools +fill-tags -Oz -o ${norm_vcf}
	vt decompose ${norm_vcf} | bcftools +fill-tags -Oz -o ${decom_vcf}
	python ~/GitCode/MitoImputePrep/scripts/PYTHON/pickFirstAlt ${decom_vcf} | bcftools view -Oz -o ${vcf_1kg}
	bcftools index ${vcf_1kg}
	plink --vcf ${vcf_1kg} --recode --double-id --keep-allele-order --out ${plink_1kg}
	bcftools query -l ${vcf_1kg} > ${samps_1kg}
	Rscript ~/GitCode/MitoImputePrep/scripts/R/assign_sex_label.R ${samps_1kg} ${sex_1kg}
	cp ${vcf_1kg} ~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/
	cp ${sex_1kg} ~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/
fi

# CREATE DIRECTORY
if [ -d  /g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/ ]
then
	echo
	echo "/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/ EXISTS ... PASSING"
else
	echo
	echo "/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/ NOT FOUND ... CREATING DIRECTORY"
	mkdir /g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/
fi

# EXTRACT PLATFORM SNPs
echo
echo "EXTRACTING PLATFORM SNPs"
#MtPlatforms=BDCHP-1X10-HUMANHAP240S_11216501_A-b37
MTSnps=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${MtPlatforms}_MT_snps.txt
vcf_1kg=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.vcf.gz
vcf=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.vcf.gz

bcftools view -R ${MTSnps} ${vcf_1kg} -Oz -o ${vcf}
bcftools index ${vcf}

# GENERATE GEN SAMPLE
echo
echo "GENERATING GEN SAMPLE"
#vcf=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}.vcf.gz
sex=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/SampleList1kg_sex.txt 
out=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}

bcftools convert --gensample ${out} ${vcf} --sex ${sex}
Rscript ~/GitCode/MitoImputePrep/scripts/R/FixSamplesFile_raijin.R ${out}.samples

# GENERATE PLINK FILES
echo
echo "GENERATING PLINK FILES"
out=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}

plink1.9 --vcf ${vcf} --recode --double-id --keep-allele-order --out ${out}

# RUN IMPUTE2
echo
echo "RUNNING IMPUTE2 ON ${MtPlatforms}"
m=~/GitCode/MitoImputePrep/DerivedData/${REFpanel}/${REFpanel}_MtMap.txt 
h=~/GitCode/MitoImputePrep/DerivedData/${REFpanel}/${REFpanel}.hap.gz
l=~/GitCode/MitoImputePrep/DerivedData/${REFpanel}/${REFpanel}.legend.gz
g=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.gen.gz
s=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.samples
out=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed

if [ -f ${out} ]
then
	echo "${out} FOUND! ... PASSING"
else
	echo "${out} NOT FOUND! ... RUNNING IMPUTE2"
	impute2 -chrX -m ${m} -h ${h} -l ${l} -g ${g} -sample_g ${s} -int 1 16569 -Ne 20000 -o ${out}
fi

# FIX CHROMOSOME NAMES
echo
echo "FIXING CHROMOSOME NAMES"
InFile=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed
OutFile=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed_ChromFixed
awk '{{$1 = "26"; print}}' ${InFile} > ${OutFile}

# CONVERT OXFORD TO PEDIGREE
echo
echo "CONVERTING OXFORD TO PEDIGREE"
gen=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed_ChromFixed
sam=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed_samples
out=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed

plink1.9 --gen ${gen} --sample ${sam} --hard-call-threshold 0.49 --keep-allele-order --output-chr 26 --recode --out ${out}

# CONVERT OXFORD TO VCF
echo
echo "CONVERTING OXFORD TO VCF"
gen=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed_ChromFixed
sam=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed_samples
out=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed

plink1.9 --gen ${gen} --sample ${sam} --hard-call-threshold 0.49 --keep-allele-order --output-chr 26 --recode vcf --out ${out}


# GENERATE QC REPORT
echo
echo "GENERATING QC REPORT"
s=~/GitCode/MitoImputePrep/scripts/R/MT_imputation_QC.Rmd
wgs_map=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.map
wgs_ped=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.ped
wgs_vcf=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.vcf.gz
typ_map=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.map
typ_ped=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.ped
typ_vcf=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.vcf.gz
imp_map=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed.map
imp_ped=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed.ped
imp_vcf=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed.vcf
imp_info=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed_info

output=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_mtImputed_QC.html

rwd = `pwd`/
output_dir=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/${MtPlatforms}/
info_cut='0'