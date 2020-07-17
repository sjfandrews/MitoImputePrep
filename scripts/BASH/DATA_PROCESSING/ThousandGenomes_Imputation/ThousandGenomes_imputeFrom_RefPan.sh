#!/bin/bash
#PBS -P gw26
#PBS -q biodev
#PBS -l walltime=00:05:00
#PBS -l mem=24GB
#PBS -l ncpus=1
#PBS -m e
#PBS -M u5015730@anu.edu.au
#PBS -j oe
#PBS -o /g/data1a/te53/MitoImpute/logs/
#PBS -l storage=scratch/te53+gdata/te53

# LOAD THE MODULE
#module unload intel-fc intel-cc
#module load python3/3.7.4
#module load python2/2.7.17
module load Rpackages/3.6.1
module load bcftools/1.9
module load plink/1.90
#module load impute2/2.3.2
module load vt
#module load java/jdk1.8.0_60
#echo
#echo "LOADED R v3.4.3"
#echo "LOADED bcftools v1.8"
#echo "LOADED plink v1.9"
#echo "LOADED IMPUTE2 v2.3.2"

echo
echo "YOU REQUIRE R v3.6.2"
echo "YOU REQUIRE bcftools v1.8"
echo "YOU REQUIRE plink v1.9"
echo "YOU REQUIRE IMPUTE2 v2.3.2"
echo

# TURN ON RESOURCE MONITORING
echo
job_num=`basename ${PBS_JOBFS} .gadi-pbs`
resource_log=/g/data1a/te53/MitoImpute/logs/${job_num}_resourceLog.txt
echo "SAVING RESOURCE LOG TO:	${resource_log}"
sh ~/GitCode/haploco/scripts/BASH/Other/monitor_job_resource_usage.sh ${job_num} ${resource_log} 300  &
first_status=`${nqstat_anu} ${job_num} | sed '2q;d'`

echo
echo "test start"
echo "FIRST STATUS UPDATE:"
echo ${first_status}
echo "test end"
echo

# ASSIGN VARIABLES FROM COMMAND LINE
REFpanel=${REFpanel}
MtPlatforms=${MtPlatforms}
mcmc=${mcmc}
burn=${burn}
khap=${khap}
ne=${ne}

if [[ "${REFpanel}" == "ReferencePanel_v2" ]] ||  [[ "${REFpanel}" == "ReferencePanel_v1-unique_0.01" ]] || [[ "${REFpanel}" == "ReferencePanel_v1_0.01" ]]
then
	maf_str="MAF1%"
	echo "MINOR ALLELE FREQUENCY = ${maf_str}!"
elif [[ "${REFpanel}" == "ReferencePanel_v4" ]] ||  [[ "${REFpanel}" == "ReferencePanel_v1-unique_0.05" ]] || [[ "${REFpanel}" == "ReferencePanel_v1_0.05" ]]
then
	maf_str="MAF0.5%"
	echo "MINOR ALLELE FREQUENCY = ${maf_str}!"
elif [[ "${REFpanel}" == "ReferencePanel_v4" ]] ||  [[ "${REFpanel}" == "ReferencePanel_v1-unique_0.05" ]] || [[ "${REFpanel}" == "ReferencePanel_v1_0.05" ]]
then
	maf_str="MAF0.1%"
	echo "MINOR ALLELE FREQUENCY = ${maf_str}!"
else
	echo "SOMETHING IS WRONG	...	EXITING!"
	exit	
fi


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



# DEFINE MAIN DIRECTORY
mitoimpute_dir=/g/data1a/te53/MitoImpute/
# DEFINE SUB-DIRECTORIES
mitoimpute_vcf_dir=${mitoimpute_dir}data/VCF/
mitoimpute_fasta_dir=${mitoimpute_dir}data/FASTA/
mitoimpute_plink_dir=${mitoimpute_dir}data/PLINK/
mitoimpute_oxford_dir=${mitoimpute_dir}data/OXFORD/
strand_dir=${mitoimpute_dir}data/STRANDS/${MtPlatforms}/${REFpanel}/
imp_dir=${strand_dir}MCMC${mcmc}/
REF_PAN_HG_DIR=~/GitCode/MitoImputePrep/metadata/HaploGrep_concordance/${REFpanel}/
thousand_g_dir=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/
REF_PAN_HG_DIR=~/GitCode/MitoImputePrep/metadata/HaploGrep_concordance/${REFpanel}/



# CHECK FOR OR CREATE THE DECOMPOSED 1,000 GENOMES VCF FILE
git_vcf=${thousand_g_dir}chrMT_1kg_norm_decomposed_firstAlt.vcf.gz
git_vcf_simple=${thousand_g_dir}chrMT_1kg_norm_decomposed_firstAlt.vcf.gz
norm_vcf=${mitoimpute_vcf_dir}chrMT_1kg_norm.vcf.gz
decom_vcf=${mitoimpute_vcf_dir}chrMT_1kg_norm_decomposed.vcf.gz
vcf_1kg=${mitoimpute_vcf_dir}chrMT_1kg_norm_decomposed_firstAlt.vcf.gz
plink_1kg=${mitoimpute_plink_dir}chrMT_1kg_norm_decomposed_firstAlt
orig_vcf=${mitoimpute_vcf_dir}ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz
snpOnly_vcf=${mitoimpute_vcf_dir}chrMT_1kg_SNPonly
diploid_vcf=${mitoimpute_vcf_dir}chrMT_1kg_diploid
samps_1kg=${mitoimpute_dir}metadata/SampleList1kg.txt
sex_1kg=${mitoimpute_dir}metadata/SampleList1kg_sex.txt
ref_fasta=${mitoimpute_fasta_dir}rCRS.fasta


if [ -f ${git_vcf} ] || [ -f ${diploid_vcf}.vcf.gz ]
then
	echo
	echo "${git_vcf}  or ${diploid_vcf}.vcf.gz EXISTS ... PASSING"
else
	echo
	echo "${git_vcf} NOT FOUND ... CREATE THE DECOMPOSED 1,000 GENOMES VCF FILE"
	
	if [ ! -d ${thousand_g_dir} ]
	then
		echo
		echo "${thousand_g_dir} NOT FOUND ... CREATING DIRECTORY"
		mkdir -p ${thousand_g_dir}
	fi
	
	bcftools norm -f ${ref_fasta} -m - ${orig_vcf} | bcftools view -V indels,mnps | bcftools norm -m + | bcftools +fill-tags -Oz -o ${norm_vcf}
	vt decompose ${norm_vcf} | bcftools +fill-tags -Oz -o ${decom_vcf}
	python3.8 ~/GitCode/MitoImputePrep/scripts/PYTHON/pickFirstAlt.py ${decom_vcf} | bcftools view -Oz -o ${vcf_1kg}
	bcftools index ${vcf_1kg}
	plink --vcf ${vcf_1kg} --recode --double-id --keep-allele-order --out ${plink_1kg}
	bcftools query -l ${vcf_1kg} > ${samps_1kg}
	Rscript ~/GitCode/MitoImputePrep/scripts/R/DATA_PROCESSING/assign_sex_label.R ${samps_1kg} ${sex_1kg}
	
	#if [ ! -f ${diploid_vcf}.vcf.gz ]
	if [ ! -f ${thousand_g_dir}chrMT_1kg_diploid.vcf.gz ]
	then
		bcftools view -V indels,mnps ${orig_vcf} | bcftools norm -m -any -Oz -o ${snpOnly_vcf}.vcf.gz
		bcftools index ${snpOnly_vcf}.vcf.gz
		plink --vcf ${snpOnly_vcf}.vcf.gz --recode vcf --out ${diploid_vcf}
		bcftools +fill-tags ${diploid_vcf}.vcf -Oz -o ${diploid_vcf}.vcf.gz
		bcftools index ${diploid_vcf}.vcf.gz
		java -jar ${HAPLOGREP} --in ${diploid_vcf}.vcf --format vcf --out ${diploid_vcf}.haplogrep.txt
		
		cp ${vcf_1kg} ${thousand_g_dir}
		cp ${vcf_1kg}.csi ${thousand_g_dir}
		cp ${sex_1kg} ${thousand_g_dir}
		cp ${diploid_vcf}.vcf.gz* ${thousand_g_dir}
		cp ${diploid_vcf}.haplogrep.txt ~/GitCode/MitoImputePrep/metadata/HaploGrep_concordance/
		
		echo
		echo "${diploid_vcf}.vcf.gz COPIED TO ${thousand_g_dir}"
	else
		echo
		echo "${diploid_vcf}.vcf.gz FOUND ... PASSING"
	fi
	
	
fi


# CREATE DIRECTORY
if [ ! -d  ${strand_dir} ]
then
	echo
	echo "${strand_dir} NOT FOUND ... CREATING DIRECTORY"
	mkdir ${strand_dir}
fi


if [ ! -d ${REF_PAN_HG_DIR} ]
then
	echo
	echo "${REF_PAN_HG_DIR} NOT FOUND ... CREATING DIRECTORY"
	mkdir -p ${REF_PAN_HG_DIR}
fi

# EXTRACT PLATFORM SNPs
#MtPlatforms=BDCHP-1X10-HUMANHAP240S_11216501_A-b37
MTSnps=${mitoimpute_dir}data/STRANDS/${MtPlatforms}/${MtPlatforms}_MT_snps.txt
vcf_1kg=${thousand_g_dir}chrMT_1kg_norm_decomposed_firstAlt.vcf.gz
vcf=${strand_dir}chrMT_1kg_${MtPlatforms}.vcf.gz

if [ ! -s ${vcf} ]
then
	echo
	echo "EXTRACTING PLATFORM SNPs"
	bcftools view -R ${MTSnps} ${vcf_1kg} -Oz -o ${vcf}
	bcftools index ${vcf}
fi


# GENERATE GEN SAMPLE

#vcf=${thousand_g_dir}${MtPlatforms}/chrMT_1kg_${MtPlatforms}.vcf.gz
sex=${thousand_g_dir}SampleList1kg_sex.txt 
out=${strand_dir}chrMT_1kg_${MtPlatforms}

if [ ! -s ${out}.gen.gz ] && [ ! -s ${out}.gen.gz ]
then
	echo
	echo "GENERATING GEN SAMPLE"
	#bcftools convert --haplegendsample ${out} --haploid2diploid ${vcf} --sex ${sex}
	bcftools convert --gensample ${out} ${vcf} --sex ${sex}
	Rscript ~/GitCode/MitoImputePrep/scripts/R/DATA_PROCESSING/FixSamplesFile_raijin.R ${out}.samples
else
	echo 
	echo "GEN AND SAMPLE FILES FOUND	...	PASSING"
fi

# GENERATE plink FILES
out=${strand_dir}chrMT_1kg_${MtPlatforms}

if [ ! -s ${out}.ped ] && [ ! -s ${out}.map ]
then
	echo
	echo "GENERATING plink FILES"
	plink --vcf ${vcf} --recode --double-id --keep-allele-order --out ${out}
else
	echo 
	echo "PED AND MAP FILES FOUND	...	PASSING"
fi


# CREATE DIPLOID VCF
geno_ext="${strand_dir}chrMT_1kg_${MtPlatforms}"
ref_fasta_plink=~/GitCode/MitoImputePrep/scripts/REFERENCE_ALNS/MT/rCRS.fasta
norm_vcf=${geno_ext}_norm.vcf.gz
vcf_pos=${geno_ext}_norm_SNPpositions.txt
geno_fasta=${geno_ext}.fasta
fixed_vcf=${geno_ext}_fixed.vcf
diploid_vcf=${geno_ext}_diploid
haplogrep_file=${diploid_vcf}_haplogrep.txt

if [ ! -s ${diploid_vcf}.vcf.gz ] && [ ! -s ${norm_vcf} ] && [ ! -s ${vcf_pos} ] && [ ! -s ${geno_fasta} ] && [ ! -s ${fixed_vcf} ]
then
	echo
	echo "GENERATING PLINK FILES (DIPLOID)"
	bcftools annotate --set-id '.' ${vcf} | bcftools norm --check-ref s -f ${ref_fasta_plink} -m +any | bcftools view -Oz -o ${norm_vcf} # Normalise: remove SNP IDs, reformat to rCRS, join biallelic repeated sites into multiallelic sites, then output to gzip
	bcftools index ${norm_vcf} # index normalised vcf
	bcftools query -f '%POS\n' ${norm_vcf} > ${vcf_pos} # extract genomic positions
	Rscript ~/GitCode/MitoImputePrep/scripts/R/DATA_PROCESSING/plink_sites_map.R ${vcf_pos} # add a column with the MT label
	perl -pi -e 'chomp if eof' ${vcf_pos} # remove the last leading line
	python2 ~/GitCode/MitoImputePrep/scripts/PYTHON/vcf2fasta_rCRS.py -i ${norm_vcf} -o ${geno_fasta} -v # convert to a fasta file
	#python ~/GitCode/MitoImputePrep/scripts/PYTHON/fasta2vcf_mtDNA.py -i ${imp_fasta} -o ${fixed_vcf} -g -d # convert back to a vcf
	python2 ~/GitCode/MitoImputePrep/scripts/PYTHON/fasta2vcf_mtDNA.py -i ${geno_fasta} -o ${fixed_vcf} -g -d -id -a -v # convert back to a vcf
	bcftools view ${fixed_vcf} -Oz -o ${fixed_vcf}.gz # gzip it so the -R flag in bcftools view will work
	bcftools index ${fixed_vcf}.gz # index it it so the -R flag in bcftools view will work
	#bcftools view -R ${vcf_pos} ${fixed_vcf}.gz | bcftools norm -m -any -Oz -o ${final_vcf}.vcf.gz # include only positions found in the imputed vcf and split multiallelic into biallelic
	bcftools view -R ${vcf_pos} ${fixed_vcf}.gz | bcftools norm -m -any | bcftools +fill-tags -Oz -o ${diploid_vcf}.vcf.gz # include only positions found in the imputed vcf and split multiallelic into biallelic
	bcftools index ${diploid_vcf}.vcf.gz # index it
else
	echo
	echo "${diploid_vcf}.vcf.gz FOUND	...	PASSING"
fi

if [ ! -s ${diploid_vcf}_haplogrep.txt ]
then
	java -jar ${HAPLOGREP} --in ${diploid_vcf}.vcf.gz --format vcf --chip --out ${haplogrep_file} # assign haplogreps
	cp ${haplogrep_file} ${REF_PAN_HG_DIR}
fi

if [ -s ${diploid_vcf}_haplogrep.txt ]
then
	echo
	echo "${diploid_vcf}.txt FOUND ... bcftools VCF FILE WORKED"
else
	echo
	echo "${haplogrep_file} NOT FOUND ... RECODING TO plink VCF FILE"
	plink --vcf ${diploid_vcf}.vcf.gz --recode vcf --out ${diploid_vcf} # recode vcf to vcf via plink (haplogrep seems to love plink vcf files, but not bcftools ... dont know why this needs to be done, but it does, so ???)
	java -jar ${HAPLOGREP} --in ${diploid_vcf}.vcf --format vcf --chip --out ${haplogrep_file} # assign haplogreps
	cp ${haplogrep_file} ${REF_PAN_HG_DIR}
fi


#plink --vcf ${vcf} --recode vcf --out ${out}
#bcftools +fill-tags ${out}.vcf -Oz -o ${out}.vcf.gz


# ASSIGN HAPLOGROUPS FOR THE PRE-IMPUTED IN SILICO SNP CHIP
echo
echo "ASSIGNING HAPLOGROUPS FOR THE PRE-IMPUTED IN SILICO SNP CHIP"
vcf="${strand_dir}chrMT_1kg_${MtPlatforms}_diploid.vcf"
out="${strand_dir}chrMT_1kg_${MtPlatforms}_haplogrep"
#java -jar ${HAPLOGREP} --in ${vcf} --format vcf --chip --out ${out}.txt # assign haplogreps
#
#if [ ! -d ~/GitCode/MitoImputePrep/metadata/HaploGrep_concordance/${REFpanel}/ ]
#then
#	mkdir -p ~/GitCode/MitoImputePrep/metadata/HaploGrep_concordance/${REFpanel}/
#fi

#cp ${out}.txt ~/GitCode/MitoImputePrep/metadata/HaploGrep_concordance/${REFpanel}/ # copy haplogroup outputs to Git

if [ -s ${haplogrep_file} ]
then
	echo "${haplogrep_file} SUCCESSFULLY COPIED"
else
	echo "HAPLOGREP FILE DID NOT COPY ... SOMETHING WENT WRONG"
fi

# RUN IMPUTE2
echo
echo "RUNNING IMPUTE2 ON ${MtPlatforms}"


m=~/GitCode/MitoImputePrep/DerivedData/${REFpanel}/${REFpanel}_MtMap.txt 
h=~/GitCode/MitoImputePrep/DerivedData/${REFpanel}/${REFpanel}.hap.gz
l=~/GitCode/MitoImputePrep/DerivedData/${REFpanel}/${REFpanel}.legend.gz
g=${strand_dir}chrMT_1kg_${MtPlatforms}.gen.gz
s=${strand_dir}chrMT_1kg_${MtPlatforms}.samples
out=${imp_dir}chrMT_1kg_${MtPlatforms}_imputed_MCMC${mcmc}

if [ -d ${imp_dir} ]
then
	echo "DIRECTORY FOUND"
else
	mkdir -p ${imp_dir}
fi

if [ -s ${out} ]
then
	echo "${out} FOUND! ... PASSING"
else
	echo "${out} NOT FOUND! ... RUNNING IMPUTE2"
	impute2 -chrX -m ${m} -h ${h} -l ${l} -g ${g} -sample_g ${s} -int 1 16569 -Ne ${ne} -o ${out} -iter ${mcmc} -burnin ${burn} -k_hap ${khap}
fi


# FIX CHROMOSOME NAMES

InFile=${imp_dir}chrMT_1kg_${MtPlatforms}_imputed_MCMC${mcmc}
OutFile=${InFile}_ChromFixed

if [ ! -s ${OutFile} ]
then
	echo
echo "FIXING CHROMOSOME NAMES"
	awk '{{$1 = "26"; print}}' ${InFile} > ${OutFile}
fi


# CONVERT OXFORD TO PEDIGREE

out_prefix=${imp_dir}chrMT_1kg_${MtPlatforms}_imputed_MCMC${mcmc}
gen=${out_prefix}_ChromFixed
sam=${out_prefix}_samples
out=${out_prefix}

if [ ! -s ${out_prefix}.ped ] && [ ! -s ${out_prefix}.map ]
then
	echo
	echo "CONVERTING OXFORD TO PEDIGREE"
	plink --gen ${gen} --sample ${sam} --hard-call-threshold 0.49 --keep-allele-order --output-chr 26 --recode --out ${out}
else
	echo 
	echo "PED AND MAP FILES FOUND	...	PASSING"
fi

# CONVERT OXFORD TO VCF
echo
echo "CONVERTING OXFORD TO VCF"
gen=${out_prefix}_ChromFixed
sam=${out_prefix}_samples
out=${out_prefix}

if [ ! -s ${out_prefix}.vcf ] 
then
	echo
	echo "CONVERTING OXFORD TO VCF"
	plink --gen ${gen} --sample ${sam} --hard-call-threshold 0.49 --keep-allele-order --output-chr 26 --recode vcf --out ${out}

else
	echo 
	echo "VCF FILES FOUND	...	PASSING"
fi


# CONVERT VCF TO FORMAT FOR HAPLOGREP2
ref_fasta_plink=~/GitCode/MitoImputePrep/scripts/REFERENCE_ALNS/26/rCRS.fasta
imp_ext=${imp_dir}chrMT_1kg_${MtPlatforms}_imputed_MCMC${mcmc}
imp_vcf=${imp_ext}.vcf
norm_imp_vcf=${imp_ext}_norm.vcf.gz
imp_fasta=${imp_ext}.fasta
vcf_pos=${imp_ext}_norm_SNPpositions.txt
fixed_vcf=${imp_ext}_fixed.vcf
final_vcf=${imp_ext}_haplogrep

if [ ! -s ${final_vcf}.txt ] && [ ! -s ${imp_vcf}.vcf ] && [ ! -s ${norm_imp_vcf} ] && [ ! -s ${imp_fasta} ] && [ ! -s ${vcf_pos} ] && [ ! -s ${fixed_vcf} ] && [ ! -s ${final_vcf}.vcf.gz ]
then
	echo
	echo "FIXING IMPUTED VCF	...	BRINGING UP TO BCFTOOLS STANDARDS AND ANNOTATING"
	echo "NORMALISED VCF SAVING TO:	${norm_imp_vcf}"
	bcftools annotate --set-id '.' ${imp_vcf} | bcftools norm --check-ref s -f ${ref_fasta_plink} -m +any | bcftools view -Oz -o ${norm_imp_vcf} # Normalise: remove SNP IDs, reformat to rCRS, join biallelic repeated sites into multiallelic sites, then output to gzip
	bcftools index ${norm_imp_vcf} # index normalised vcf
	
	echo "VCF POSITION INFO SAVING TO:	${vcf_pos}" 
	bcftools query -f '%POS\n' ${norm_imp_vcf} > ${vcf_pos} # extract genomic positions
	Rscript ~/GitCode/MitoImputePrep/scripts/R/DATA_PROCESSING/plink_sites_map.R ${vcf_pos} # add a column with the MT label
	perl -pi -e 'chomp if eof' ${vcf_pos} # remove the last leading line
	
	echo "FASTA FILE SAVING TO:	${imp_fasta}" 
	python2 ~/GitCode/MitoImputePrep/scripts/PYTHON/vcf2fasta_rCRS.py -i ${norm_imp_vcf} -o ${imp_fasta} # convert to a fasta file
	#python ~/GitCode/MitoImputePrep/scripts/PYTHON/fasta2vcf_mtDNA.py -i ${imp_fasta} -o ${fixed_vcf} -g -d # convert back to a vcf
	echo "FIXED VCF SAVING TO:	${fixed_vcf}" 
	python2 ~/GitCode/MitoImputePrep/scripts/PYTHON/fasta2vcf_mtDNA.py -i ${imp_fasta} -o ${fixed_vcf} -g -d -id -a # convert back to a vcf
	bcftools view ${fixed_vcf} -Oz -o ${fixed_vcf}.gz # gzip it so the -R flag in bcftools view will work
	bcftools index ${fixed_vcf}.gz # index it it so the -R flag in bcftools view will work
	#bcftools view -R ${vcf_pos} ${fixed_vcf}.gz | bcftools norm -m -any -Oz -o ${final_vcf}.vcf.gz # include only positions found in the imputed vcf and split multiallelic into biallelic
	echo "FINAL ANNOTATED HAPLOGREP COMPATIBLE VCF SAVING TO:	${final_vcf}.vcf.gz"
	bcftools view -R ${vcf_pos} ${fixed_vcf}.gz | bcftools norm -m -any | bcftools +fill-tags -Oz -o ${final_vcf}.vcf.gz # include only positions found in the imputed vcf and split multiallelic into biallelic
	bcftools index ${final_vcf}.vcf.gz # index it
	echo "HAPLOGREP FILE SAVING TO:	${final_vcf}.txt" 
	java -jar ${HAPLOGREP} --in ${final_vcf}.vcf.gz --format vcf --chip --out ${final_vcf}.txt # assign haplogreps
else
	echo
	echo "THE FOLLOWING FILE WERE SUPPOSEDLY FOUND:"
	echo ${final_vcf}.txt
	echo ${imp_vcf}.vcf
	echo ${norm_imp_vcf}
	echo ${imp_fasta}
	echo ${vcf_pos}
	echo ${fixed_vcf}
	echo ${final_vcf}.vcf.gz
fi

if [ -s ${final_vcf}.txt ]
then
	echo
	echo "${final_vcf}.txt FOUND ... bcftools VCF FILE WORKED"
else
	echo
	echo "${final_vcf}.txt NOT FOUND ... RECODING TO plink VCF FILE"
	plink --vcf ${final_vcf}.vcf.gz --recode vcf --out ${final_vcf} # recode vcf to vcf via plink (haplogrep seems to love plink vcf files, but not bcftools ... dont know why this needs to be done, but it does, so ???)
	java -jar ${HAPLOGREP} --in ${final_vcf}.vcf --format vcf --chip --out ${final_vcf}.txt # assign haplogreps
fi

if [ -s ${final_vcf}.txt ]
then
	echo
	echo "${final_vcf}.txt FOUND ... GREAT!"
else
	echo
	echo "${final_vcf}.txt NOT FOUND ... SOMETHING HAS GONE WRONG"
fi


## CALCULATE Matthew's Correlation Coefficient
REF26=~/GitCode/MitoImputePrep/scripts/REFERENCE_ALNS/26/rCRS.fasta 
REFMT=~/GitCode/MitoImputePrep/scripts/REFERENCE_ALNS/MT/rCRS.fasta 
WGS_VCF=${vcf_1kg}
TYP_VCF=${strand_dir}chrMT_1kg_${MtPlatforms}.vcf.gz
TYP_VCF_MULTIALLELIC=${strand_dir}chrMT_1kg_${MtPlatforms}_multiallelic.vcf.gz
TYP_VCF_DECOMPOSED=${strand_dir}chrMT_1kg_${MtPlatforms}_biallelic_decomposed.vcf.gz
#IMP_VCF=${imp_ext}.vcf
IMP_INFO=${imp_ext}_info
OUT_FILE=${imp_ext}

if [ ! -s ${TYP_VCF_DECOMPOSED} ]
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

mcc_imputed=${OUT_FILE}_imputed_MCC.csv
mcc_typed=${OUT_FILE}_typed_MCC.csv

if [ ! -s ${IMP_VCF_rCRS} ]
then
	bcftools norm --check-ref s -f ${REF26} -m + ${imp_ext}.vcf -Oz -o ${IMP_VCF_rCRS}
fi

# DECOMPOSE IMPUTED VCF
if [ ! -s ${IMP_VCF} ]
then
	vt decompose ${IMP_VCF_rCRS} | bcftools +fill-tags -Oz -o ${IMP_VCF}
	bcftools index ${IMP_VCF}
fi

if [ -s ${mcc_imputed} ] & [ -f ${mcc_typed} ]
then
	echo
	echo "${mcc_imputed} AND ${mcc_typed} FOUND ... PIPELINE COMPLETED"
	echo "ACTUALLY ... DO IT ANYWAY (DOUBLE CHECKING, REMOVE THIS LATER"
	#Rscript ~/GitCode/MitoImputePrep/scripts/R/ANALYSIS/MCC/MCC_Genotypes.R ${WGS_VCF} ${TYP_VCF_DECOMPOSED} ${IMP_VCF} ${IMP_INFO} ${OUT_FILE}
else
	echo
	echo "${mcc_imputed} AND ${mcc_typed} NOT FOUND ... CALCULATING MCC GENOTYPE CONCORDANCE"
	Rscript ~/GitCode/MitoImputePrep/scripts/R/ANALYSIS/MCC/MCC_Genotypes.R ${WGS_VCF} ${TYP_VCF_DECOMPOSED} ${IMP_VCF} ${IMP_INFO} ${OUT_FILE}
fi


# CUTOFF BY IMPUTE2 INFO SCORE
impute2_file=${imp_dir}chrMT_1kg_${MtPlatforms}_imputed_MCMC${mcmc}
impute2_info_file=${impute2_file}_info
impute2_file_cutoff=${impute2_file}_cutoffRetained
impute2_info_file_cutoff=${impute2_info_file}_cutoffRetained
cutoff=0.3

if [ ! -s ${impute2_file_cutoff} ] && [ ! -s ${impute2_info_file_cutoff} ]
then
	echo
	echo "${impute2_file_cutoff} AND ${impute2_info_file_cutoff} NOT FOUND	...	CUTTING OFF AT INFO ≥ ${cutoff}"
	Rscript ~/GitCode/MitoImputePrep/scripts/R/DATA_PROCESSING/remove_impute2_cutoff.R ${impute2_file}
else
	echo
	echo "${impute2_file_cutoff} AND ${impute2_info_file_cutoff} FOUND	...	PASSING"
fi


# FIX CHROMOSOME NAMES

InFile=${impute2_file_cutoff}
OutFile=${InFile}_ChromFixed

if [ ! -s ${OutFile} ]
then
	echo
	echo "FIXING CHROMOSOME NAMES"
	awk '{{$1 = "26"; print}}' ${InFile} > ${OutFile}
fi


# CONVERT OXFORD TO PEDIGREE

out_prefix_cutoff=${impute2_file_cutoff}
gen=${out_prefix_cutoff}_ChromFixed
sam=${out_prefix}_samples
out=${out_prefix_cutoff}

if [ ! -s ${out_prefix_cutoff}.ped ] && [ ! -s ${out_prefix_cutoff}.map ]
then
	echo
	echo "CONVERTING OXFORD TO PEDIGREE"
	plink --gen ${gen} --sample ${sam} --hard-call-threshold 0.49 --keep-allele-order --output-chr 26 --recode --out ${out}
else
	echo 
	echo "PED AND MAP FILES FOUND	...	PASSING"
fi


# CONVERT OXFORD TO VCF
echo
echo "CONVERTING OXFORD TO VCF"
gen=${out_prefix_cutoff}_ChromFixed
sam=${out_prefix}_samples
out=${out_prefix_cutoff}

if [ ! -s ${out_prefix_cutoff}.vcf ] 
then
	echo
	echo "CONVERTING OXFORD TO VCF"
	plink --gen ${gen} --sample ${sam} --hard-call-threshold 0.49 --keep-allele-order --output-chr 26 --recode vcf --out ${out}

else
	echo 
	echo "VCF FILES FOUND	...	PASSING"
fi

# CONVERT VCF TO FORMAT FOR HAPLOGREP2
ref_fasta_plink=~/GitCode/MitoImputePrep/scripts/REFERENCE_ALNS/26/rCRS.fasta
imp_ext=${out_prefix_cutoff}
imp_vcf=${imp_ext}.vcf
norm_imp_vcf=${imp_ext}_norm.vcf.gz
imp_fasta=${imp_ext}.fasta
vcf_pos=${imp_ext}_norm_SNPpositions.txt
fixed_vcf=${imp_ext}_fixed.vcf
final_vcf=${imp_ext}_haplogrep

if [ -s ${imp_vcf}.vcf ] & [ -s ${norm_imp_vcf} ] & [ -s ${imp_fasta} ] & [ -s ${vcf_pos} ] & [ -s ${fixed_vcf} ] & [ -s ${final_vcf}.vcf.gz ]
then
	echo
	echo "THE FOLLOWING FILE WERE SUPPOSEDLY FOUND:"
	echo ${final_vcf}.txt
	echo ${imp_vcf}.vcf
	echo ${norm_imp_vcf}
	echo ${imp_fasta}
	echo ${vcf_pos}
	echo ${fixed_vcf}
	echo ${final_vcf}.vcf.gz
else
	echo
	echo "FIXING IMPUTED VCF	...	BRINGING UP TO BCFTOOLS STANDARDS AND ANNOTATING"
	echo "NORMALISED VCF SAVING TO:	${norm_imp_vcf}"
	bcftools annotate --set-id '.' ${imp_vcf} | bcftools norm --check-ref s -f ${ref_fasta_plink} -m +any | bcftools view -Oz -o ${norm_imp_vcf} # Normalise: remove SNP IDs, reformat to rCRS, join biallelic repeated sites into multiallelic sites, then output to gzip
	bcftools index ${norm_imp_vcf} # index normalised vcf
	
	echo "VCF POSITION INFO SAVING TO:	${vcf_pos}" 
	bcftools query -f '%POS\n' ${norm_imp_vcf} > ${vcf_pos} # extract genomic positions
	Rscript ~/GitCode/MitoImputePrep/scripts/R/DATA_PROCESSING/plink_sites_map.R ${vcf_pos} # add a column with the MT label
	perl -pi -e 'chomp if eof' ${vcf_pos} # remove the last leading line
	
	echo "FASTA FILE SAVING TO:	${imp_fasta}" 
	python2 ~/GitCode/MitoImputePrep/scripts/PYTHON/vcf2fasta_rCRS.py -i ${norm_imp_vcf} -o ${imp_fasta} # convert to a fasta file
	#python ~/GitCode/MitoImputePrep/scripts/PYTHON/fasta2vcf_mtDNA.py -i ${imp_fasta} -o ${fixed_vcf} -g -d # convert back to a vcf
	echo "FIXED VCF SAVING TO:	${fixed_vcf}" 
	python2 ~/GitCode/MitoImputePrep/scripts/PYTHON/fasta2vcf_mtDNA.py -i ${imp_fasta} -o ${fixed_vcf} -g -d -id -a # convert back to a vcf
	bcftools view ${fixed_vcf} -Oz -o ${fixed_vcf}.gz # gzip it so the -R flag in bcftools view will work
	bcftools index ${fixed_vcf}.gz # index it it so the -R flag in bcftools view will work
	#bcftools view -R ${vcf_pos} ${fixed_vcf}.gz | bcftools norm -m -any -Oz -o ${final_vcf}.vcf.gz # include only positions found in the imputed vcf and split multiallelic into biallelic
	echo "FINAL ANNOTATED HAPLOGREP COMPATIBLE VCF SAVING TO:	${final_vcf}.vcf.gz"
	bcftools view -R ${vcf_pos} ${fixed_vcf}.gz | bcftools norm -m -any | bcftools +fill-tags -Oz -o ${final_vcf}.vcf.gz # include only positions found in the imputed vcf and split multiallelic into biallelic
	bcftools index ${final_vcf}.vcf.gz # index it
	#echo "HAPLOGREP FILE SAVING TO:	${final_vcf}.txt" 
	#java -jar ${HAPLOGREP} --in ${final_vcf}.vcf.gz --format vcf --chip --out ${final_vcf}.txt # assign haplogreps
fi

if [ -s ${final_vcf}.txt ]
then
	echo
	echo "${final_vcf}.txt FOUND ... bcftools VCF FILE WORKED"
else
	echo "HAPLOGREP FILE SAVING TO:	${final_vcf}.txt" 
	java -jar ${HAPLOGREP} --in ${final_vcf}.vcf.gz --format vcf --chip --out ${final_vcf}.txt # assign haplogreps
fi


if [ -s ${final_vcf}.txt ]
then
	echo
	echo "${final_vcf}.txt FOUND ... bcftools VCF FILE WORKED"
else
	echo
	echo "${final_vcf}.txt NOT FOUND ... RECODING TO plink VCF FILE"
	plink --vcf ${final_vcf}.vcf.gz --recode vcf --out ${final_vcf} # recode vcf to vcf via plink (haplogrep seems to love plink vcf files, but not bcftools ... dont know why this needs to be done, but it does, so ???)
	java -jar ${HAPLOGREP} --in ${final_vcf}.vcf --format vcf --chip --out ${final_vcf}.txt # assign haplogreps
fi

exit
exit
exit

## CALCULATE Matthew's Correlation Coefficient
REF26=~/GitCode/MitoImputePrep/scripts/REFERENCE_ALNS/26/rCRS.fasta 
REFMT=~/GitCode/MitoImputePrep/scripts/REFERENCE_ALNS/MT/rCRS.fasta 
WGS_VCF=${vcf_1kg}
TYP_VCF=${strand_dir}chrMT_1kg_${MtPlatforms}.vcf.gz
TYP_VCF_MULTIALLELIC=${strand_dir}chrMT_1kg_${MtPlatforms}_multiallelic.vcf.gz
TYP_VCF_DECOMPOSED=${strand_dir}chrMT_1kg_${MtPlatforms}_biallelic_decomposed.vcf.gz
#IMP_VCF=${imp_ext}.vcf
IMP_INFO=${impute2_info_file_cutoff}
OUT_FILE=${out_prefix_cutoff}
IMP_VCF_rCRS=${out_prefix_cutoff}_rCRS.vcf
IMP_VCF=${out_prefix_cutoff}_FINAL.vcf.gz

if [ ! -s ${TYP_VCF_DECOMPOSED} ]
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


if [ ! -s ${IMP_VCF_rCRS} ]
then
	bcftools norm --check-ref s -f ${REF26} -m + ${imp_ext}.vcf -Oz -o ${IMP_VCF_rCRS}
fi

# DECOMPOSE IMPUTED VCF
if [ ! -s ${IMP_VCF} ]
then
	vt decompose ${IMP_VCF_rCRS} | bcftools +fill-tags -Oz -o ${IMP_VCF}
	bcftools index ${IMP_VCF}
fi

mcc_imputed_cutoff=${OUT_FILE}_imputed_MCC.csv
mcc_typed_cutoff=${OUT_FILE}_typed_MCC.csv

if [ -s ${mcc_imputed_cutoff} ] & [ -f ${mcc_typed_cutoff} ]
then
	echo
	echo "${mcc_imputed_cutoff} AND ${mcc_typed_cutoff} FOUND ... PIPELINE COMPLETED"
	echo "ACTUALLY ... DO IT ANYWAY (DOUBLE CHECKING, REMOVE THIS LATER"
	#Rscript ~/GitCode/MitoImputePrep/scripts/R/ANALYSIS/MCC/MCC_Genotypes.R ${WGS_VCF} ${TYP_VCF_DECOMPOSED} ${IMP_VCF} ${IMP_INFO} ${OUT_FILE}
else
	echo
	echo "${mcc_imputed_cutoff} AND ${mcc_typed_cutoff} NOT FOUND ... CALCULATING MCC GENOTYPE CONCORDANCE"
	Rscript ~/GitCode/MitoImputePrep/scripts/R/ANALYSIS/MCC/MCC_Genotypes.R ${WGS_VCF} ${TYP_VCF_DECOMPOSED} ${IMP_VCF} ${IMP_INFO} ${OUT_FILE}
fi

# GENERATE HiMC HAPLOGROUPINGS
imp_ext=${imp_dir}chrMT_1kg_${MtPlatforms}_imputed_MCMC${mcmc}
full_1kGP_pref=${mitoimpute_dir}data/PLINK/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes
typed_1kGP_pref=${strand_dir}chrMT_1kg_${MtPlatforms}
imputed_1kGP_pref=${imp_ext}
imputed_cutoff_1kGP_pref=${imp_ext}_cutoffRetained

himc_hg_file=${imp_ext}_HiMC_haplogroups.csv


if [ -s ${himc_hg_file} ]
then
	echo
	echo "${himc_hg_file} FOUND!	...	PASSING"
else
	echo
	echo "${himc_hg_file} NOT FOUND!	...	GENERATING HiMC HAPLOGROUPINGS"
	Rscript ~/GitCode/MitoImputePrep/scripts/R/DATA_PROCESSING/HiMC_haplogrouping.R ${full_1kGP_pref} ${typed_1kGP_pref} ${imputed_1kGP_pref} ${imputed_cutoff_1kGP_pref}
fi


# GENERATE HAPLOGREP HAPLOGROUPINGS
imp_ext=${imp_dir}chrMT_1kg_${MtPlatforms}_imputed_MCMC${mcmc}
full_1kGP_hg=~/GitCode/MitoImputePrep/DerivedData/ThousandGenomes/1000genomes_mtDNA_haplogrep.txt
typed_1kGP_hg=${strand_dir}chrMT_1kg_${MtPlatforms}_diploid_haplogrep.txt
imputed_1kGP_hg=${imp_ext}_haplogrep.txt
imputed_cutoff_1kGP_hg=${imp_ext}_cutoffRetained_haplogrep.txt

haplogrep_hg_file=${imp_ext}_HaploGrep_haplogroups.csv


if [ -s ${haplogrep_hg_file} ]
then
	echo
	echo "${haplogrep_hg_file} FOUND!	...	PASSING"
else
	echo
	echo "${haplogrep_hg_file} NOT FOUND!	...	GENERATING HiMC HAPLOGROUPINGS"
	Rscript ~/GitCode/MitoImputePrep/scripts/R/DATA_PROCESSING/HaploGrep_haplogrouping.R ${full_1kGP_hg} ${typed_1kGP_hg} ${imputed_1kGP_hg} ${imputed_cutoff_1kGP_hg}
fi


# SUMMARISE EVERYTHING!
final_summary_file=${imp_dir}chrMT_1kg_${MtPlatforms}_imputed_MCMC${mcmc}_SUMMARY.csv

mcmc_str="MCMC${mcmc}"
khap_str="kHAP${khap}"

ls -lh ${final_summary_file}
echo
echo ${MtPlatforms}
echo ${mcmc_str}
echo ${maf_str}
echo ${khap_str}
echo ${himc_hg_file}
echo ${haplogrep_hg_file}
echo ${impute2_info_file}
echo ${impute2_info_file_cutoff}
echo ${MTSnps}
echo ${mcc_imputed}
echo ${mcc_imputed_cutoff}
echo


if [ -s ${final_summary_file} ]
then
	echo
	echo "${final_summary_file} FOUND	...	PASSING"
else
	echo
	echo "${final_summary_file} NOT FOUND	...	SUMMARISING EVERYTHING"
	Rscript ~/GitCode/MitoImputePrep/scripts/R/DATA_PROCESSING/summarise_HiMC_HaploGrep_perchip.R ${MtPlatforms} ${mcmc_str} ${maf_str} ${khap_str} ${himc_hg_file} ${haplogrep_hg_file} ${impute2_info_file} ${impute2_info_file_cutoff} ${MTSnps} ${mcc_imputed} ${mcc_imputed_cutoff}
fi 

echo
echo "END!"
echo

# END!