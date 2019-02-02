#!/bin/bash
#PBS -P te53
#PBS -q normalbw
#PBS -l walltime=48:00:00
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
REFpanel="ReferencePanel_v2"
HAPLOGREP=~/GitCode/MitoImputePrep/haplogrep/2.1.19/haplogrep-2.1.19.jar
echo
echo "REFERENCE PANEL:	${REFpanel}"
echo "SNP CHIP:			${MtPlatforms}"
echo "HAPLOGREP:		${HAPLOGREP}"


# CHECK FOR OR CREATE THE DECOMPOSED 1,000 GENOMES VCF FILE
git_vcf=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.vcf.gz
norm_vcf=/g/data1a/te53/MitoImpute/data/VCF/chrMT_1kg_norm.vcf.gz
decom_vcf=/g/data1a/te53/MitoImpute/data/VCF/chrMT_1kg_norm_decomposed.vcf.gz
vcf_1kg=/g/data1a/te53/MitoImpute/data/VCF/chrMT_1kg_norm_decomposed_firstAlt.vcf.gz
plink_1kg=/g/data1a/te53/MitoImpute/data/PLINK/chrMT_1kg_norm_decomposed_firstAlt
orig_vcf=/g/data1a/te53/haploco/data/originals/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz
samps_1kg=/g/data1a/te53/MitoImpute/metadata/SampleList1kg.txt
sex_1kg=/g/data1a/te53/MitoImpute/metadata/SampleList1kg_sex.txt
ref_fasta=/g/data1a/te53/MitoImpute/data/FASTA/rCRS.fasta
if [ -f ${git_vcf} ]
then
	echo
	echo "${git_vcf} EXISTS ... PASSING"
else
	echo
	echo "${git_vcf} NOT FOUND ... CREATE THE DECOMPOSED 1,000 GENOMES VCF FILE"
	bcftools norm -f ${ref_fasta} -m - ${orig_vcf} | bcftools view -V indels,mnps | bcftools norm -m + | bcftools +fill-tags -Oz -o ${norm_vcf}
	vt decompose ${norm_vcf} | bcftools +fill-tags -Oz -o ${decom_vcf}
	python ~/GitCode/MitoImputePrep/scripts/PYTHON/pickFirstAlt ${decom_vcf} | bcftools view -Oz -o ${vcf_1kg}
	bcftools index ${vcf_1kg}
	plink --vcf ${vcf_1kg} --recode --double-id --keep-allele-order --out ${plink_1kg}
	bcftools query -l ${vcf_1kg} > ${samps_1kg}
	Rscript ~/GitCode/MitoImputePrep/scripts/R/assign_sex_label.R ${samps_1kg} ${sex_1kg}
	cp ${vcf_1kg} ~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/
	cp ${vcf_1kg}.csi ~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/
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

# CREATE DIPLOID VCF
echo
echo "GENERATING PLINK FILES (DIPLOID)"
geno_ext="/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}"
ref_fasta_plink=~/GitCode/MitoImputePrep/scripts/REFERENCE_ALNS/MT/rCRS.fasta
norm_vcf=${geno_ext}_norm.vcf.gz
vcf_pos=${geno_ext}_norm_SNPpositions.txt
geno_fasta=${geno_ext}.fasta
fixed_vcf=${geno_ext}_fixed.vcf
diploid_vcf=${geno_ext}_diploid

bcftools annotate --set-id '.' ${vcf} | bcftools norm --check-ref s -f ${ref_fasta_plink} -m +any | bcftools view -Oz -o ${norm_vcf} # Normalise: remove SNP IDs, reformat to rCRS, join biallelic repeated sites into multiallelic sites, then output to gzip
bcftools index ${norm_vcf} # index normalised vcf
bcftools query -f '%POS\n' ${norm_vcf} > ${vcf_pos} # extract genomic positions
Rscript ~/GitCode/MitoImputePrep/scripts/R/plink_sites_map.R ${vcf_pos} # add a column with the MT label
perl -pi -e 'chomp if eof' ${vcf_pos} # remove the last leading line
python ~/GitCode/MitoImputePrep/scripts/PYTHON/vcf2fasta_rCRS.py -i ${norm_vcf} -o ${geno_fasta} # convert to a fasta file
#python ~/GitCode/MitoImputePrep/scripts/PYTHON/fasta2vcf_mtDNA.py -i ${imp_fasta} -o ${fixed_vcf} -g -d # convert back to a vcf
python ~/GitCode/MitoImputePrep/scripts/PYTHON/fasta2vcf_mtDNA.py -i ${geno_fasta} -o ${fixed_vcf} -g -d -id -a # convert back to a vcf
bcftools view ${fixed_vcf} -Oz -o ${fixed_vcf}.gz # gzip it so the -R flag in bcftools view will work
bcftools index ${fixed_vcf}.gz # index it it so the -R flag in bcftools view will work
#bcftools view -R ${vcf_pos} ${fixed_vcf}.gz | bcftools norm -m -any -Oz -o ${final_vcf}.vcf.gz # include only positions found in the imputed vcf and split multiallelic into biallelic
bcftools view -R ${vcf_pos} ${fixed_vcf}.gz | bcftools norm -m -any | bcftools +fill-tags -Oz -o ${diploid_vcf}.vcf.gz # include only positions found in the imputed vcf and split multiallelic into biallelic
bcftools index ${diploid_vcf}.vcf.gz # index it

java -jar ${HAPLOGREP} --in ${diploid_vcf}.vcf.gz --format vcf --chip --out ${diploid_vcf}.txt # assign haplogreps

if [ -f ${diploid_vcf}.txt ]
then
	echo
	echo "${diploid_vcf}.txt FOUND ... bcftools VCF FILE WORKED"
else
	echo
	echo "${diploid_vcf}.txt NOT FOUND ... RECODING TO plink VCF FILE"
	plink1.9 --vcf ${diploid_vcf}.vcf.gz --recode vcf --out ${diploid_vcf} # recode vcf to vcf via plink (haplogrep seems to love plink vcf files, but not bcftools ... dont know why this needs to be done, but it does, so ???)
	java -jar ${HAPLOGREP} --in ${diploid_vcf}.vcf --format vcf --chip --out ${diploid_vcf}.txt # assign haplogreps
fi

# RUN IMPUTE2
echo
echo "RUNNING IMPUTE2 ON ${MtPlatforms}"
m=~/GitCode/MitoImputePrep/DerivedData/${REFpanel}/${REFpanel}_MtMap.txt 
h=~/GitCode/MitoImputePrep/DerivedData/${REFpanel}/${REFpanel}.hap.gz
l=~/GitCode/MitoImputePrep/DerivedData/${REFpanel}/${REFpanel}.legend.gz
g=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.gen.gz
s=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.samples
out=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/kHAP_Experiments/kHAP${khap}/chrMT_1kg_${MtPlatforms}_imputed_kHAP${khap}

if [ -d /g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/kHAP_Experiments/kHAP${khap}/ ]
then
	echo "DIRECTORY FOUND"
else
	mkdir -p /g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/kHAP_Experiments/kHAP${khap}/
fi

if [ -f ${out} ]
then
	echo "${out} FOUND! ... PASSING"
else
	echo "${out} NOT FOUND! ... RUNNING IMPUTE2"
	impute2 -chrX -m ${m} -h ${h} -l ${l} -g ${g} -sample_g ${s} -int 1 16569 -Ne 20000 -o ${out} -iter ${mcmc} -burnin ${burn} -k_hap ${khap}
fi

# FIX CHROMOSOME NAMES
echo
echo "FIXING CHROMOSOME NAMES"
InFile=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/kHAP_Experiments/kHAP${khap}/chrMT_1kg_${MtPlatforms}_imputed_kHAP${khap}
OutFile=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/kHAP_Experiments/kHAP${khap}/chrMT_1kg_${MtPlatforms}_imputed_kHAP${khap}_ChromFixed
awk '{{$1 = "26"; print}}' ${InFile} > ${OutFile}

# CONVERT OXFORD TO PEDIGREE
echo
echo "CONVERTING OXFORD TO PEDIGREE"
gen=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/kHAP_Experiments/kHAP${khap}/chrMT_1kg_${MtPlatforms}_imputed_kHAP${khap}_ChromFixed
sam=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/kHAP_Experiments/kHAP${khap}/chrMT_1kg_${MtPlatforms}_imputed_kHAP${khap}_samples
out=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/kHAP_Experiments/kHAP${khap}/chrMT_1kg_${MtPlatforms}_imputed_kHAP${khap}

plink1.9 --gen ${gen} --sample ${sam} --hard-call-threshold 0.49 --keep-allele-order --output-chr 26 --recode --out ${out}

# CONVERT OXFORD TO VCF
echo
echo "CONVERTING OXFORD TO VCF"
gen=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/kHAP_Experiments/kHAP${khap}/chrMT_1kg_${MtPlatforms}_imputed_kHAP${khap}_ChromFixed
sam=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/kHAP_Experiments/kHAP${khap}/chrMT_1kg_${MtPlatforms}_imputed_kHAP${khap}_samples
out=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/kHAP_Experiments/kHAP${khap}/chrMT_1kg_${MtPlatforms}_imputed_kHAP${khap}

plink1.9 --gen ${gen} --sample ${sam} --hard-call-threshold 0.49 --keep-allele-order --output-chr 26 --recode vcf --out ${out}

# CONVERT VCF TO FORMAT FOR HAPLOGREP2
ref_fasta_plink=~/GitCode/MitoImputePrep/scripts/REFERENCE_ALNS/26/rCRS.fasta
imp_ext=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/kHAP_Experiments/kHAP${khap}/chrMT_1kg_${MtPlatforms}_imputed_kHAP${khap}
imp_vcf=${imp_ext}.vcf
norm_imp_vcf=${imp_ext}_norm.vcf.gz
imp_fasta=${imp_ext}.fasta
vcf_pos=${imp_ext}_norm_SNPpositions.txt
fixed_vcf=${imp_ext}_fixed.vcf
final_vcf=${imp_ext}_haplogrep

bcftools annotate --set-id '.' ${imp_vcf} | bcftools norm --check-ref s -f ${ref_fasta_plink} -m +any | bcftools view -Oz -o ${norm_imp_vcf} # Normalise: remove SNP IDs, reformat to rCRS, join biallelic repeated sites into multiallelic sites, then output to gzip
bcftools index ${norm_imp_vcf} # index normalised vcf
bcftools query -f '%POS\n' ${norm_imp_vcf} > ${vcf_pos} # extract genomic positions
Rscript ~/GitCode/MitoImputePrep/scripts/R/plink_sites_map.R ${vcf_pos} # add a column with the MT label
perl -pi -e 'chomp if eof' ${vcf_pos} # remove the last leading line
python ~/GitCode/MitoImputePrep/scripts/PYTHON/vcf2fasta_rCRS.py -i ${norm_imp_vcf} -o ${imp_fasta} # convert to a fasta file
#python ~/GitCode/MitoImputePrep/scripts/PYTHON/fasta2vcf_mtDNA.py -i ${imp_fasta} -o ${fixed_vcf} -g -d # convert back to a vcf
python ~/GitCode/MitoImputePrep/scripts/PYTHON/fasta2vcf_mtDNA.py -i ${imp_fasta} -o ${fixed_vcf} -g -d -id -a # convert back to a vcf
bcftools view ${fixed_vcf} -Oz -o ${fixed_vcf}.gz # gzip it so the -R flag in bcftools view will work
bcftools index ${fixed_vcf}.gz # index it it so the -R flag in bcftools view will work
#bcftools view -R ${vcf_pos} ${fixed_vcf}.gz | bcftools norm -m -any -Oz -o ${final_vcf}.vcf.gz # include only positions found in the imputed vcf and split multiallelic into biallelic
bcftools view -R ${vcf_pos} ${fixed_vcf}.gz | bcftools norm -m -any | bcftools +fill-tags -Oz -o ${final_vcf}.vcf.gz # include only positions found in the imputed vcf and split multiallelic into biallelic
bcftools index ${final_vcf}.vcf.gz # index it
java -jar ${HAPLOGREP} --in ${final_vcf}.vcf.gz --format vcf --chip --out ${final_vcf}.txt # assign haplogreps

if [ -f ${final_vcf}.txt ]
then
	echo
	echo "${final_vcf}.txt FOUND ... bcftools VCF FILE WORKED"
else
	echo
	echo "${final_vcf}.txt NOT FOUND ... RECODING TO plink VCF FILE"
	plink1.9 --vcf ${final_vcf}.vcf.gz --recode vcf --out ${final_vcf} # recode vcf to vcf via plink (haplogrep seems to love plink vcf files, but not bcftools ... dont know why this needs to be done, but it does, so ???)
	java -jar ${HAPLOGREP} --in ${final_vcf}.vcf --format vcf --chip --out ${final_vcf}.txt # assign haplogreps
fi


# GENERATE QC REPORT
#echo
#echo "GENERATING QC REPORT"
#s=~/GitCode/MitoImputePrep/scripts/R/MT_imputation_QC.Rmd
#wgs_map=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.map
#wgs_ped=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.ped
#wgs_vcf=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.vcf.gz
#typ_map=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.map
#typ_ped=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.ped
#typ_vcf=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.vcf.gz
#imp_map=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed.map
#imp_ped=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed.ped
#imp_vcf=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed.vcf
#imp_info=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed_info
#
#output=/g/data1a/te53/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_mtImputed_QC.html
#
#rwd = `pwd`/
#output_dir=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/${MtPlatforms}/
#info_cut='0'