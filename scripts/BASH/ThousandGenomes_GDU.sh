#!/bin/bash
#$ -m e
#$ -M u5015730@anu.edu.au
#$ -b y
#$ -l h_vmem=5g,virtual_free=4.9g
#$ -N impute_SNPchip_1kGP
#$ -o /home/easteallab/tim/MitoImpute/logs/
#$ -j y
#set -ex

# SPECIFY REFERENCE PANEL
REFpanel="ReferencePanel_v3"
echo
echo "REFERENCE PANEL:	${REFpanel}"
echo "SNP CHIP:			${MtPlatforms}"

# CHECK FOR OR CREATE THE DECOMPOSED 1,000 GENOMES VCF FILE
git_vcf=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.vcf.gz
norm_vcf=/home/easteallab/tim/MitoImpute/data/VCF/chrMT_1kg_norm.vcf.gz
decom_vcf=/home/easteallab/tim/MitoImpute/data/VCF/chrMT_1kg_norm_decomposed.vcf.gz
vcf_1kg=/home/easteallab/tim/MitoImpute/data/VCF/chrMT_1kg_norm_decomposed_firstAlt.vcf.gz
plink_1kg=/home/easteallab/tim/MitoImpute/data/PLINK/chrMT_1kg_norm_decomposed_firstAlt
orig_vcf=/home/easteallab/tim/haploco/data/originals/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz
samps_1kg=/home/easteallab/tim/MitoImpute/metadata/SampleList1kg.txt
sex_1kg=/home/easteallab/tim/MitoImpute/metadata/SampleList1kg_sex.txt
ref_fasta=/home/easteallab/tim/MitoImpute/data/FASTA/rCRS.fasta
if [ -f ${git_vcf} ]
then
	echo
	echo "${git_vcf} EXISTS ... PASSING"
else
	echo
	echo "${git_vcf} NOT FOUND ... CREATE THE DECOMPOSED 1,000 GENOMES VCF FILE"
	/home/easteallab/software/bcftools norm -f ${ref_fasta} -m - ${orig_vcf} | /home/easteallab/software/bcftools view -V indels,mnps | /home/easteallab/software/bcftools norm -m + | /home/easteallab/software/bcftools +fill-tags -Oz -o ${norm_vcf}
	vt decompose ${norm_vcf} | /home/easteallab/software/bcftools +fill-tags -Oz -o ${decom_vcf}
	python ~/GitCode/MitoImputePrep/scripts/PYTHON/pickFirstAlt ${decom_vcf} | /home/easteallab/software/bcftools view -Oz -o ${vcf_1kg}
	/home/easteallab/software/bcftools index ${vcf_1kg}
	plink1.9 --vcf ${vcf_1kg} --recode --double-id --keep-allele-order --out ${plink_1kg}
	/home/easteallab/software/bcftools query -l ${vcf_1kg} > ${samps_1kg}
	Rscript ~/GitCode/MitoImputePrep/scripts/R/assign_sex_label.R ${samps_1kg} ${sex_1kg}
	cp ${vcf_1kg} ~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/
	cp ${vcf_1kg}.csi ~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/
	cp ${sex_1kg} ~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/
fi

# CREATE DIRECTORY
if [ -d  /home/easteallab/tim/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/ ]
then
	echo
	echo "/home/easteallab/tim/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/ EXISTS ... PASSING"
else
	echo
	echo "/home/easteallab/tim/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/ NOT FOUND ... CREATING DIRECTORY"
	mkdir /home/easteallab/tim/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/
fi

# EXTRACT PLATFORM SNPs
echo
echo "EXTRACTING PLATFORM SNPs"
#MtPlatforms=BDCHP-1X10-HUMANHAP240S_11216501_A-b37
MTSnps=/home/easteallab/tim/MitoImpute/data/STRANDS/${MtPlatforms}/${MtPlatforms}_MT_snps.txt
vcf_1kg=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.vcf.gz
vcf=/home/easteallab/tim/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.vcf.gz

/home/easteallab/software/bcftools view -R ${MTSnps} ${vcf_1kg} -Oz -o ${vcf}
/home/easteallab/software/bcftools index ${vcf}

# GENERATE GEN SAMPLE
echo
echo "GENERATING GEN SAMPLE"
#vcf=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}.vcf.gz
sex=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/SampleList1kg_sex.txt 
out=/home/easteallab/tim/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}

/home/easteallab/software/bcftools convert --gensample ${out} ${vcf} --sex ${sex}
Rscript ~/GitCode/MitoImputePrep/scripts/R/FixSamplesFile_raijin.R ${out}.samples

# GENERATE PLINK FILES
echo
echo "GENERATING PLINK FILES"
out=/home/easteallab/tim/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}

/home/easteallab/software/plink1.9 --vcf ${vcf} --recode --double-id --keep-allele-order --out ${out}

# RUN IMPUTE2
echo
echo "RUNNING IMPUTE2 ON ${MtPlatforms}"
m=~/GitCode/MitoImputePrep/DerivedData/${REFpanel}/${REFpanel}_MtMap.txt 
h=~/GitCode/MitoImputePrep/DerivedData/${REFpanel}/${REFpanel}.hap.gz
l=~/GitCode/MitoImputePrep/DerivedData/${REFpanel}/${REFpanel}.legend.gz
g=/home/easteallab/tim/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.gen.gz
s=/home/easteallab/tim/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.samples
out=/home/easteallab/tim/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed

if [ -f ${out} ]
then
	echo "${out} FOUND! ... PASSING"
else
	echo "${out} NOT FOUND! ... RUNNING IMPUTE2"
	/home/easteallab/software/impute2 -chrX -m ${m} -h ${h} -l ${l} -g ${g} -sample_g ${s} -int 1 16569 -Ne 20000 -o ${out}
fi

# FIX CHROMOSOME NAMES
echo
echo "FIXING CHROMOSOME NAMES"
InFile=/home/easteallab/tim/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed
OutFile=/home/easteallab/tim/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed_ChromFixed
awk '{{$1 = "26"; print}}' ${InFile} > ${OutFile}

# CONVERT OXFORD TO PEDIGREE
echo
echo "CONVERTING OXFORD TO PEDIGREE"
gen=/home/easteallab/tim/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed_ChromFixed
sam=/home/easteallab/tim/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed_samples
out=/home/easteallab/tim/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed

/home/easteallab/software/plink1.9 --gen ${gen} --sample ${sam} --hard-call-threshold 0.49 --keep-allele-order --output-chr 26 --recode --out ${out}

# CONVERT OXFORD TO VCF
echo
echo "CONVERTING OXFORD TO VCF"
gen=/home/easteallab/tim/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed_ChromFixed
sam=/home/easteallab/tim/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed_samples
out=/home/easteallab/tim/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed

/home/easteallab/software/plink1.9 --gen ${gen} --sample ${sam} --hard-call-threshold 0.49 --keep-allele-order --output-chr 26 --recode vcf --out ${out}

echo
echo "######################################################################################################"
echo "											PROCESS COMPLETE											"
echo "######################################################################################################"
echo