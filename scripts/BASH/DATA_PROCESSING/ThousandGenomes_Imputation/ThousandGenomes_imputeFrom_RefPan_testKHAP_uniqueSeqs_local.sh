#!/bin/bash


echo
echo "YOU REQUIRE R v3.6.2"
echo "YOU REQUIRE bcftools v1.8"
echo "YOU REQUIRE plink1.9 v1.9"
echo "YOU REQUIRE IMPUTE2 v2.3.2"

# ASSIGN VARIABLES FROM COMMAND LINE
REFpanel=$1
MtPlatforms=$2
mcmc=$3
burn=$4
khap=$5
ne=$6

# SPECIFY REFERENCE PANEL
#REFpanel="ReferencePanel_v5"
REFpanel=${REFpanel}
HAPLOGREP=~/GitCode/MitoImputePrep/haplogrep/2.1.25/haplogrep-2.1.25.jar
echo
echo "REFERENCE PANEL:				${REFpanel}"
echo "SNP CHIP:					${MtPlatforms}"
echo "HAPLOGREP:					${HAPLOGREP}"
echo "MCMC LENGTH:					${mcmc}"
echo "MCMC BURN-IN:					${burn}"
echo "REF HAPLOTYPES (k_hap):				${khap}"
echo "EFFECTIVE POP. SIZE (Ne):			${ne}"
echo

# CHECK FOR OR CREATE THE DECOMPOSED 1,000 GENOMES VCF FILE
git_vcf=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.vcf.gz
git_vcf_simple=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.vcf.gz
norm_vcf=/Volumes/TimMcInerney/MitoImpute/data/VCF/chrMT_1kg_norm.vcf.gz
decom_vcf=/Volumes/TimMcInerney/MitoImpute/data/VCF/chrMT_1kg_norm_decomposed.vcf.gz
vcf_1kg=/Volumes/TimMcInerney/MitoImpute/data/VCF/chrMT_1kg_norm_decomposed_firstAlt.vcf.gz
plink_1kg=/Volumes/TimMcInerney/MitoImpute/data/PLINK/chrMT_1kg_norm_decomposed_firstAlt
orig_vcf=/Volumes/TimMcInerney/MitoImpute/data/VCF/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz
snpOnly_vcf=/Volumes/TimMcInerney/MitoImpute/data/VCF/chrMT_1kg_SNPonly
diploid_vcf=/Volumes/TimMcInerney/MitoImpute/data/VCF/chrMT_1kg_diploid
samps_1kg=/Volumes/TimMcInerney/MitoImpute/metadata/SampleList1kg.txt
sex_1kg=/Volumes/TimMcInerney/MitoImpute/metadata/SampleList1kg_sex.txt
ref_fasta=/Volumes/TimMcInerney/MitoImpute/data/FASTA/rCRS.fasta

if [ -f ${git_vcf} ] || [ -f ${diploid_vcf}.vcf.gz ]
then
	echo
	echo "${git_vcf}  or ${diploid_vcf}.vcf.gz EXISTS ... PASSING"
else
	echo
	echo "${git_vcf} NOT FOUND ... CREATE THE DECOMPOSED 1,000 GENOMES VCF FILE"
	
	if [ ! -d ~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/ ]
	then
		echo
		echo "~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/ NOT FOUND ... CREATING DIRECTORY"
		mkdir -p ~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/
	fi
	
	bcftools norm -f ${ref_fasta} -m - ${orig_vcf} | bcftools view -V indels,mnps | bcftools norm -m + | bcftools +fill-tags -Oz -o ${norm_vcf}
	vt decompose ${norm_vcf} | bcftools +fill-tags -Oz -o ${decom_vcf}
	python3.8 ~/GitCode/MitoImputePrep/scripts/PYTHON/pickFirstAlt.py ${decom_vcf} | bcftools view -Oz -o ${vcf_1kg}
	bcftools index ${vcf_1kg}
	plink1.9 --vcf ${vcf_1kg} --recode --double-id --keep-allele-order --out ${plink_1kg}
	bcftools query -l ${vcf_1kg} > ${samps_1kg}
	Rscript ~/GitCode/MitoImputePrep/scripts/R/DATA_PROCESSING/assign_sex_label.R ${samps_1kg} ${sex_1kg}
	
	#if [ ! -f ${diploid_vcf}.vcf.gz ]
	if [ ! -f ~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/chrMT_1kg_diploid.vcf.gz ]
	then
		bcftools view -V indels,mnps ${orig_vcf} | bcftools norm -m -any -Oz -o ${snpOnly_vcf}.vcf.gz
		bcftools index ${snpOnly_vcf}.vcf.gz
		plink1.9 --vcf ${snpOnly_vcf}.vcf.gz --recode vcf --out ${diploid_vcf}
		bcftools +fill-tags ${diploid_vcf}.vcf -Oz -o ${diploid_vcf}.vcf.gz
		bcftools index ${diploid_vcf}.vcf.gz
		java -jar ${HAPLOGREP} --in ${diploid_vcf}.vcf --format vcf --out ${diploid_vcf}.haplogrep.txt
		
		cp ${vcf_1kg} ~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/
		cp ${vcf_1kg}.csi ~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/
		cp ${sex_1kg} ~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/
		cp ${diploid_vcf}.vcf.gz* ~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/
		cp ${diploid_vcf}.haplogrep.txt ~/GitCode/MitoImputePrep/metadata/HaploGrep_concordance/
		
		echo
		echo "${diploid_vcf}.vcf.gz COPIED TO ~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/"
	else
		echo
		echo "${diploid_vcf}.vcf.gz FOUND ... PASSING"
	fi
	
	
fi


# CREATE DIRECTORY
if [ -d  /Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/ ]
then
	echo
	echo "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/ EXISTS ... PASSING"
else
	echo
	echo "/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/ NOT FOUND ... CREATING DIRECTORY"
	mkdir /Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/
fi

# EXTRACT PLATFORM SNPs
echo
echo "EXTRACTING PLATFORM SNPs"
#MtPlatforms=BDCHP-1X10-HUMANHAP240S_11216501_A-b37
MTSnps=/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${MtPlatforms}_MT_snps.txt
vcf_1kg=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.vcf.gz
vcf=/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.vcf.gz

bcftools view -R ${MTSnps} ${vcf_1kg} -Oz -o ${vcf}
bcftools index ${vcf}
exit
# GENERATE GEN SAMPLE
echo
echo "GENERATING GEN SAMPLE"
#vcf=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/${MtPlatforms}/chrMT_1kg_${MtPlatforms}.vcf.gz
sex=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/SampleList1kg_sex.txt 
out=/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}

#bcftools convert --haplegendsample ${out} --haploid2diploid ${vcf} --sex ${sex}
bcftools convert --gensample ${out} ${vcf} --sex ${sex}
Rscript ~/GitCode/MitoImputePrep/scripts/R/FixSamplesFile_raijin.R ${out}.samples

# GENERATE plink1.9 FILES
echo
echo "GENERATING plink1.9 FILES"
out=/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}

plink1.9 --vcf ${vcf} --recode --double-id --keep-allele-order --out ${out}

# CREATE DIPLOID VCF
echo
echo "GENERATING plink1.9 FILES (DIPLOID)"
geno_ext="/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}"
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
	echo "${diploid_vcf}.txt NOT FOUND ... RECODING TO plink1.9 VCF FILE"
	plink1.9 --vcf ${diploid_vcf}.vcf.gz --recode vcf --out ${diploid_vcf} # recode vcf to vcf via plink1.9 (haplogrep seems to love plink1.9 vcf files, but not bcftools ... dont know why this needs to be done, but it does, so ???)
	java -jar ${HAPLOGREP} --in ${diploid_vcf}.vcf --format vcf --chip --out ${diploid_vcf}.txt # assign haplogreps
fi

#plink1.9 --vcf ${vcf} --recode vcf --out ${out}
#bcftools +fill-tags ${out}.vcf -Oz -o ${out}.vcf.gz

# ASSIGN HAPLOGROUPS FOR THE PRE-IMPUTED IN SILICO SNP CHIP
echo
echo "ASSIGNING HAPLOGROUPS FOR THE PRE-IMPUTED IN SILICO SNP CHIP"
vcf="/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_diploid.vcf"
out="/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_haplogrep"
#java -jar ${HAPLOGREP} --in ${vcf} --format vcf --chip --out ${out}.txt # assign haplogreps
#
#if [ ! -d ~/GitCode/MitoImputePrep/metadata/HaploGrep_concordance/${REFpanel}/ ]
#then
#	mkdir -p ~/GitCode/MitoImputePrep/metadata/HaploGrep_concordance/${REFpanel}/
#fi

#cp ${out}.txt ~/GitCode/MitoImputePrep/metadata/HaploGrep_concordance/${REFpanel}/ # copy haplogroup outputs to Git

if [ -f ~/GitCode/MitoImputePrep/metadata/HaploGrep_concordance/${REFpanel}/chrMT_1kg_${MtPlatforms}_haplogrep.txt ]
then
	echo "~/GitCode/MitoImputePrep/metadata/HaploGrep_concordance/${REFpanel}/chrMT_1kg_${MtPlatforms}_haplogrep.txt SUCCESSFULLY COPIED"
else
	echo "HAPLOGREP FILE DID NOT COPY ... SOMETHING WENT WRONG"
fi


# RUN IMPUTE2
echo
echo "RUNNING IMPUTE2 ON ${MtPlatforms}"
m=~/GitCode/MitoImputePrep/DerivedData/${REFpanel}/${REFpanel}_MtMap.txt 
h=~/GitCode/MitoImputePrep/DerivedData/${REFpanel}/${REFpanel}.hap.gz
l=~/GitCode/MitoImputePrep/DerivedData/${REFpanel}/${REFpanel}.legend.gz
g=/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.gen.gz
s=/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.samples
out=/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/kHAP${khap}_uniqueseqs/chrMT_1kg_${MtPlatforms}_imputed_kHAP${khap}

if [ -d /Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/kHAP${khap}_uniqueseqs/ ]
then
	echo "DIRECTORY FOUND"
else
	mkdir -p /Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/kHAP${khap}_uniqueseqs/
fi

if [ -f ${out} ]
then
	echo "${out} FOUND! ... PASSING"
else
	echo "${out} NOT FOUND! ... RUNNING IMPUTE2"
	impute2 -chrX -m ${m} -h ${h} -l ${l} -g ${g} -sample_g ${s} -int 1 16569 -Ne ${ne} -o ${out} -iter ${mcmc} -burnin ${burn} -k_hap ${khap}
fi

# FIX CHROMOSOME NAMES
echo
echo "FIXING CHROMOSOME NAMES"
InFile=/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/kHAP${khap}_uniqueseqs/chrMT_1kg_${MtPlatforms}_imputed_kHAP${khap}
OutFile=/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/kHAP${khap}_uniqueseqs/chrMT_1kg_${MtPlatforms}_imputed_kHAP${khap}_ChromFixed
awk '{{$1 = "26"; print}}' ${InFile} > ${OutFile}

# CONVERT OXFORD TO PEDIGREE
echo
echo "CONVERTING OXFORD TO PEDIGREE"
gen=/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/kHAP${khap}_uniqueseqs/chrMT_1kg_${MtPlatforms}_imputed_kHAP${khap}_ChromFixed
sam=/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/kHAP${khap}_uniqueseqs/chrMT_1kg_${MtPlatforms}_imputed_kHAP${khap}_samples
out=/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/kHAP${khap}_uniqueseqs/chrMT_1kg_${MtPlatforms}_imputed_kHAP${khap}

plink1.9 --gen ${gen} --sample ${sam} --hard-call-threshold 0.49 --keep-allele-order --output-chr 26 --recode --out ${out}

# CONVERT OXFORD TO VCF
echo
echo "CONVERTING OXFORD TO VCF"
gen=/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/kHAP${khap}_uniqueseqs/chrMT_1kg_${MtPlatforms}_imputed_kHAP${khap}_ChromFixed
sam=/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/kHAP${khap}_uniqueseqs/chrMT_1kg_${MtPlatforms}_imputed_kHAP${khap}_samples
out=/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/kHAP${khap}_uniqueseqs/chrMT_1kg_${MtPlatforms}_imputed_kHAP${khap}

plink1.9 --gen ${gen} --sample ${sam} --hard-call-threshold 0.49 --keep-allele-order --output-chr 26 --recode vcf --out ${out}

# CONVERT VCF TO FORMAT FOR HAPLOGREP2
ref_fasta_plink=~/GitCode/MitoImputePrep/scripts/REFERENCE_ALNS/26/rCRS.fasta
imp_ext=/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/kHAP${khap}_uniqueseqs/chrMT_1kg_${MtPlatforms}_imputed_kHAP${khap}
imp_vcf=${imp_ext}.vcf
norm_imp_vcf=${imp_ext}_norm.vcf.gz
imp_fasta=${imp_ext}.fasta
vcf_pos=${imp_ext}_norm_SNPpositions.txt
fixed_vcf=${imp_ext}_fixed.vcf
final_vcf=${imp_ext}_haplogrep

if [ ! -f ${final_vcf}.txt ]
then
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
fi

if [ -f ${final_vcf}.txt ]
then
	echo
	echo "${final_vcf}.txt FOUND ... bcftools VCF FILE WORKED"
else
	echo
	echo "${final_vcf}.txt NOT FOUND ... RECODING TO plink1.9 VCF FILE"
	plink1.9 --vcf ${final_vcf}.vcf.gz --recode vcf --out ${final_vcf} # recode vcf to vcf via plink1.9 (haplogrep seems to love plink1.9 vcf files, but not bcftools ... dont know why this needs to be done, but it does, so ???)
	java -jar ${HAPLOGREP} --in ${final_vcf}.vcf --format vcf --chip --out ${final_vcf}.txt # assign haplogreps
fi

#if [ ! -d ~/GitCode/MitoImputePrep/metadata/HaploGrep_concordance/${REFpanel}/kHAP${khap}_uniqueseqs/ ]
#then
#	mkdir -p ~/GitCode/MitoImputePrep/metadata/HaploGrep_concordance/${REFpanel}/kHAP${khap}_uniqueseqs/
#fi
#cp ${final_vcf}.txt ~/GitCode/MitoImputePrep/metadata/HaploGrep_concordance/${REFpanel}/kHAP${khap}_uniqueseqs/ # copy haplogroup outputs to Git

if [ -f ${final_vcf}.txt ]
then
	echo
	echo "${final_vcf}.txt FOUND ... PIPELINE COMPLETED"
else
	echo
	echo "${final_vcf}.txt NOT FOUND ... SOMETHING HAS GONE WRONG"
fi

## CALCULATE Matthew's Correlation Coefficient
REF26=~/GitCode/MitoImputePrep/scripts/REFERENCE_ALNS/26/rCRS.fasta 
REFMT=~/GitCode/MitoImputePrep/scripts/REFERENCE_ALNS/MT/rCRS.fasta 
WGS_VCF=${vcf_1kg}
TYP_VCF=/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.vcf.gz
TYP_VCF_MULTIALLELIC=/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_multiallelic.vcf.gz
TYP_VCF_DECOMPOSED=/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_biallelic_decomposed.vcf.gz
#IMP_VCF=${imp_ext}.vcf
IMP_INFO=${imp_ext}_info
OUT_FILE=${imp_ext}

if [ ! -f ${TYP_VCF_DECOMPOSED} ]
then
	echo
	echo "NORMALISE GENOTYPED VCF"
	# NORMALISE GENOTYPED VCF
	bcftools norm --check-ref s -f ${REFMT} -m + ${TYP_VCF} -Oz -o ${TYP_VCF_MULTIALLELIC}
	
	# DECOMPOSE GENOTYPED VCF
	vt decompose ${TYP_VCF_MULTIALLELIC} | bcftools +fill-tags -Oz -o ${TYP_VCF_DECOMPOSED}
	bcftools index ${TYP_VCF_DECOMPOSED}
else
	echo
	echo "NORMALISED GENOTYPED VCF FOUND ... MOVING ON"
fi

echo
echo "NORMALISE IMPUTED VCF"
# NORMALISE IMPUTED VCF
IMP_VCF_rCRS=${imp_ext}_rCRS.vcf
IMP_VCF=${imp_ext}_FINAL.vcf.gz
bcftools norm --check-ref s -f ${REF26} -m + ${imp_ext}.vcf -Oz -o ${IMP_VCF_rCRS}

# DECOMPOSE IMPUTED VCF
vt decompose ${IMP_VCF_rCRS} | bcftools +fill-tags -Oz -o ${IMP_VCF}
bcftools index ${IMP_VCF}

if [ -f ${OUT_FILE}_imputed_MCC.csv ] & [ -f ${OUT_FILE}_typed_MCC.csv ]
then
	echo
	echo "${OUT_FILE}_imputed_MCC.csv AND ${OUT_FILE}_typed_MCC.csv FOUND ... PIPELINE COMPLETED"
	echo "ACTUALLY ... DO IT ANYWAY (DOUBLE CHECKING, REMOVE THIS LATER"
	Rscript ~/GitCode/MitoImputePrep/scripts/R/MCC_Genotypes.R ${WGS_VCF} ${TYP_VCF_DECOMPOSED} ${IMP_VCF} ${IMP_INFO} ${OUT_FILE}
else
	echo
	echo "${OUT_FILE}_imputed_MCC.csv AND ${OUT_FILE}_typed_MCC.csv NOT FOUND ... CALCULATING MCC GENOTYPE CONCORDANCE"
	Rscript ~/GitCode/MitoImputePrep/scripts/R/MCC_Genotypes.R ${WGS_VCF} ${TYP_VCF_DECOMPOSED} ${IMP_VCF} ${IMP_INFO} ${OUT_FILE}
fi

# GENERATE QC REPORT
#echo
#echo "GENERATING QC REPORT"
#s=~/GitCode/MitoImputePrep/scripts/R/MT_imputation_QC.Rmd
#wgs_map=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.map
#wgs_ped=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.ped
#wgs_vcf=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/chrMT_1kg_norm_decomposed_firstAlt.vcf.gz
#typ_map=/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.map
#typ_ped=/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.ped
#typ_vcf=/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}.vcf.gz
#imp_map=/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed.map
#imp_ped=/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed.ped
#imp_vcf=/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed.vcf
#imp_info=/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_imputed_info
#
#output=/Volumes/TimMcInerney/MitoImpute/data/STRANDS/${MtPlatforms}/${REFpanel}/chrMT_1kg_${MtPlatforms}_mtImputed_QC.html
#
#rwd = `pwd`/
#output_dir=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/${MtPlatforms}/
#info_cut='0'