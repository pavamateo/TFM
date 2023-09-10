#!/bin/bash
#SBATCH -p long
#SBATCH --job-name="WES_vc"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=out.%x.%j
#SBATCH --error=err.%x.%j

### check the cpus to use. Using 12 cpus in the vhio cluster reduces the number of jobs you could run simultaneously in the server. 

### USER: this script is intented to be in the same folder where the bams and all results are gonna be. It is recommended that both normal and tumoral bam files are in the same folder.

### USER in VHIO SERVER: sbatch WES_vc_norm_vcfs_v1.0.sh   *TUMOR*.recal.bam     *NORMAL(SALIVA/PBMC)*.recal.bam 


RANDOM=$(date +%s%N | cut -b10-19)
RANDOMSTAMP='temp.'$RANDOM$RANDOM$RANDOM$RANDOM


### FUNCTIONS

check_file ()
	{	
	if [ ! -f $1 ]; then
		echo "$(date)	### ERROR: File $1 does not exist. EXITING.";  exit
	fi

	if [ ! -s $1 ]; then
		echo "$(date)	### ERROR: File $1 is empty. EXITING.";  exit
	fi
	}

	

reference='/mnt/bioinfnas/prostate/daguilar/references/GRCh38_exomes/Homo_sapiens_assembly38.fasta'
intervals='/mnt/bioinfnas/prostate/daguilar/references/GRCh38_exomes/KAPA_HyperExome_capture_targets.bed'
intervalsStrelka='/mnt/bioinfnas/prostate/daguilar/references/GRCh38_exomes/KAPA_HyperExome_capture_targets.bed.gz'
vcf='/mnt/bioinfnas/prostate/daguilar/references/gatk-best-practices/hg38/somatic-hg38_af-only-gnomad.hg38.vcf.gz'

check_file ${reference}
check_file ${intervals}
check_file ${intervalsStrelka}
check_file ${vcf}


### INPUT FILES

arg1=$1
arg2=$2

if [ -z ${arg1} ]; then
	echo "$(date)	### ERROR: the name of the tumor/exo file must be passed as argument #1. EXITING.";  exit
fi

if [ -z ${arg2} ]; then
	echo "$(date)	### ERROR: the name of the normal file must be passed as argument #2. EXITING.";  exit
fi

fileTumor=$1
fileNormal=$2


nameTumor=${fileTumor%*.recal.bam}
nameNormal=${fileNormal%*.recal.bam}

if [ -z ${nameTumor} ]; then
	echo "$(date)	### ERROR: cannot parse sample name from file ${fileTumor}. EXITING.";  exit
fi

if [ -z ${nameNormal} ]; then
	echo "$(date)	### ERROR: cannot parse sample name from file ${nameNormal}. EXITING.";  exit
fi

check_file ${fileTumor}
check_file ${fileNormal}

BaseName=${fileTumor%*.recal.bam}

echo "** File tumor: "${fileTumor}" (prefix: "${nameTumor}")"
echo "** File normal: "${fileNormal}" (prefix: "${nameNormal}")"
echo "** Basename for output files: ${BaseName}"



if [ -f ./${BaseName}.bams/${fileTumor} ]; then
	echo "$(date)	### ERROR: File ./${BaseName}.bams/${fileTumor} found. Was this sample already analyzed? EXITING.";  exit
fi

######################################################################################################################################################################################################
### STEP 1: VERIFY AND PRE-PROCESS FILES
######################################################################################################################################################################################################

if true; then

	# Changing bam header for mutect2

        date
        echo "** Changing bam header of ${fileTumor} to ${nameTumor}"
        samtools view -H ${fileTumor}  | sed "s/SM:[^\t]*/SM:${nameTumor}/g" | samtools reheader - ${fileTumor} > ${RANDOMSTAMP}.bam
        check_file ${RANDOMSTAMP}.bam
        mv ${RANDOMSTAMP}.bam ${fileTumor}

	date
        echo "** Indexing ${nameTumor}"
        java -Xmx32g -jar /mnt/bioinfnas/prostate/daguilar/soft/picard/picard.jar BuildBamIndex I=${fileTumor}
	check_file ${nameTumor}.recal.bai

	date
	echo "** Changing bam header of ${fileNormal} to ${nameNormal}"
	samtools view -H ${fileNormal}  | sed "s/SM:[^\t]*/SM:${nameNormal}/g" | samtools reheader - ${fileNormal} > ${RANDOMSTAMP}.bam
	check_file ${RANDOMSTAMP}.bam
	mv ${RANDOMSTAMP}.bam ${fileNormal}

	date
	echo "** Indexing ${fileNormal}"
	java -Xmx32g -jar /mnt/bioinfnas/prostate/daguilar/soft/picard/picard.jar BuildBamIndex I=${fileNormal}
	check_file ${nameNormal}.recal.bai

	# Run pileup due to Varscan requires a pileup file from each of the bam files. This generates the files named ${fileTumor}.intervals.pileup and ${fileNormal}.intervals.pileup

	date
	echo "** Running pileup on ${fileTumor}"
	if [ -f ${fileTumor}.intervals.pileup ]; then
		echo "** File ${fileTumor}.intervals.pileup already exists. Skipping."
	else
		samtools mpileup -C 50 -B -q 1 -Q 15 --positions ${intervals} -f ${reference} ${fileTumor} > ${fileTumor}.intervals.pileup
		check_file ${fileTumor}.intervals.pileup
	fi

	date
	echo "** Running pileup on ${fileNormal}"
	if [ -f ${fileNormal}.intervals.pileup ]; then
		echo "** File ${fileNormal}.intervals.pileup already exists. Skipping."
	else
		samtools mpileup -C 50 -B -q 1 -Q 15 --positions ${intervals} -f ${reference} ${fileNormal} > ${fileNormal}.intervals.pileup
		check_file ${fileNormal}.intervals.pileup
	fi

	# mv files to own folder

	mkdir ${BaseName}_stats

	mv -f *.stats.txt ./${BaseName}_stats/
	mv -f *.ValidateSamFile.txt ./${BaseName}_stats/
	
	date
	echo "* END VERIFY AND PRE-PROCESS FILES"
fi

######################################################################################################################################################################################################
### STEP 2: MUTECT2
######################################################################################################################################################################################################

if true; then
	echo "*"
	echo "* BEGIN MUTECT2"
	echo "*"

	date
	echo "** Running Mutect2"

	#intervals in .bed file:
	/mnt/bioinfnas/prostate/daguilar/soft/gatk-4.2.5.0/gatk --java-options '-Xmx32G' Mutect2 --reference ${reference} --input ${fileTumor} --input ${fileNormal} --tumor ${nameTumor} --normal ${nameNormal} --output ${fileTumor}.Mutect.unfiltered.vcf --germline-resource ${vcf} -L ${intervals} --dont-use-soft-clipped-bases --f1r2-tar-gz ${fileTumor}.f1r2.tar.gz --native-pair-hmm-threads 12

	#no intervals:
	#/mnt/bioinfnas/prostate/daguilar/soft/gatk-4.2.5.0/gatk --java-options '-Xmx32G' Mutect2 --reference ${reference} --input ${fileTumor} --input ${fileNormal} --tumor ${nameTumor} --normal ${nameNormal} --output ${fileTumor}.Mutect.unfiltered.vcf --germline-resource ${vcf} --dont-use-soft-clipped-bases --f1r2-tar-gz ${fileTumor}.f1r2.tar.gz --native-pair-hmm-threads 16

	check_file ${fileTumor}.Mutect.unfiltered.vcf
	check_file ${fileTumor}.f1r2.tar.gz

	# LearnReadOrientationModel (advisable particularly for FFPE samples)
	/mnt/bioinfnas/prostate/daguilar/soft/gatk-4.2.5.0/gatk --java-options '-Xmx32G' LearnReadOrientationModel -I ${fileTumor}.f1r2.tar.gz -O ${fileTumor}.read-orientation-model.tar.gz
	check_file ${fileTumor}.read-orientation-model.tar.gz

	date
	echo "** Filtering and selecting Variants for FFPE samples"

	# FilterMutectCalls
	/mnt/bioinfnas/prostate/daguilar/soft/gatk-4.2.5.0/gatk --java-options '-Xmx32G' FilterMutectCalls --variant ${fileTumor}.Mutect.unfiltered.vcf  --ob-priors ${fileTumor}.read-orientation-model.tar.gz --output ${fileTumor}.temp.vcf --reference ${reference}
	check_file ${fileTumor}.temp.vcf

	# SelectVariants
	/mnt/bioinfnas/prostate/daguilar/soft/gatk-4.2.5.0/gatk --java-options '-Xmx32G' SelectVariants -V ${fileTumor}.temp.vcf -O ${fileTumor}.Mutect.passed.vcf --exclude-filtered
	check_file ${fileTumor}.Mutect.passed.vcf

	# CountVariants
	/mnt/bioinfnas/prostate/daguilar/soft/gatk-4.2.5.0/gatk --java-options '-Xmx32G' CountVariants -V ${fileTumor}.Mutect.passed.vcf -VS LENIENT > ${fileTumor}.Mutect.passed.CountVariants.txt
	check_file ${fileTumor}.Mutect.passed.CountVariants.txt

	# remove temporary files
	rm -f ${fileTumor}.temp.vcf
	rm -f ${fileTumor}.temp.vcf.idx
	rm -f ${fileTumor}.temp.vcf.filteringStats.tsv
	rm -f ${fileTumor}.read-orientation-model.tar.gz
	rm -f ${fileTumor}.f1r2.tar.gz

	# mv files to own folder
	mkdir ${BaseName}_mutect
	mv -f *.Mutect.* ./${BaseName}_mutect/

	date
	echo "*"
	echo "* END MUTECT2"
	echo "*"
fi

######################################################################################################################################################################################################
### STEP 3: VARSCAN
######################################################################################################################################################################################################

if true; then
	echo "*"
	echo "* BEGIN VARSCAN"
	echo "*"

	date
	echo "** Running VarScan"
	# with intervals:
	java -jar /mnt/bioinfnas/prostate/daguilar/soft/VarScan.v2.4.3.jar somatic ${fileNormal}.intervals.pileup ${fileTumor}.intervals.pileup ${fileTumor}.varscan --tumor-purity .5 --output-vcf 1 --min-coverage 4 --min-var-freq .05 --min-reads 2 --strand-filter 1

	# no intervals:
	#java -jar /mnt/bioinfnas/prostate/daguilar/soft/VarScan.v2.4.3.jar somatic ${fileNormal}.pileup ${fileTumor}.pileup ${fileTumor}.varscan --tumor-purity .5 --output-vcf 1 --min-coverage 4 --min-var-freq .05 --min-reads 2 --strand-filter 1

	check_file ${fileTumor}.varscan.snp.vcf
	check_file ${fileTumor}.varscan.indel.vcf

	date
	echo "** Filtering by confidence (SNPs)"

	java -jar /mnt/bioinfnas/prostate/daguilar/soft/VarScan.v2.4.3.jar processSomatic ${fileTumor}.varscan.snp.vcf
	check_file ${fileTumor}.varscan.snp.Somatic.vcf
	check_file ${fileTumor}.varscan.snp.Somatic.hc.vcf
	check_file ${fileTumor}.varscan.snp.Germline.vcf
	check_file ${fileTumor}.varscan.snp.Germline.hc.vcf

	date
	echo "** Filtering by confidence (INDELs)"

	java -jar /mnt/bioinfnas/prostate/daguilar/soft/VarScan.v2.4.3.jar processSomatic ${fileTumor}.varscan.indel.vcf
	check_file ${fileTumor}.varscan.indel.Somatic.vcf
	check_file ${fileTumor}.varscan.indel.Somatic.hc.vcf
	check_file ${fileTumor}.varscan.indel.Germline.vcf
	check_file ${fileTumor}.varscan.indel.Germline.hc.vcf

	# mv files to own folder
	mkdir ${BaseName}_varscan
	mv -f *.varscan.* ./${BaseName}_varscan/
	echo "*"
	echo "* END VARSCAN"
	echo "*"
fi

######################################################################################################################################################################################################
### STEP 4: STRELKA
######################################################################################################################################################################################################

if true; then
	echo "*"
	echo "* BEGIN STRELKA"
	echo "*"
	date
	echo "** Running Strelka"

	# no intervals:
	#/mnt/bioinfnas/prostate/daguilar/soft/Python-2.7.18/python /mnt/bioinfnas/prostate/daguilar/soft/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py --normalBam ${fileNormal} --tumorBam ${fileTumor} --referenceFasta ${reference} --exome --runDir ${fileTumor}.strelka2.temp

	# with intervals:
	/mnt/bioinfnas/prostate/daguilar/soft/Python-2.7.18/python /mnt/bioinfnas/prostate/daguilar/soft/strelka-2.9.2.centos6_x86_64/bin/configureStrelkaSomaticWorkflow.py --normalBam ${fileNormal} --tumorBam ${fileTumor} --referenceFasta ${reference} --exome --callRegions ${intervalsStrelka} --runDir ${fileTumor}.strelka2.temp

	/mnt/bioinfnas/prostate/daguilar/soft/Python-2.7.18/python ${fileTumor}.strelka2.temp/runWorkflow.py -m local -j 15

	check_file ${fileTumor}.strelka2.temp/results/variants/somatic.snvs.vcf.gz
	check_file ${fileTumor}.strelka2.temp/results/variants/somatic.indels.vcf.gz

	date
	echo "** Unzipping files"
	gunzip -kf ${fileTumor}.strelka2.temp/results/variants/somatic.snvs.vcf.gz
	gunzip -kf ${fileTumor}.strelka2.temp/results/variants/somatic.indels.vcf.gz

	check_file ${fileTumor}.strelka2.temp/results/variants/somatic.snvs.vcf
	check_file ${fileTumor}.strelka2.temp/results/variants/somatic.indels.vcf

	mv ${fileTumor}.strelka2.temp/results/variants/somatic.snvs.vcf ${fileTumor}.strelka2.snp.Somatic.vcf 
	mv ${fileTumor}.strelka2.temp/results/variants/somatic.indels.vcf ${fileTumor}.strelka2.indel.Somatic.vcf

	check_file ${fileTumor}.strelka2.snp.Somatic.vcf 
	check_file ${fileTumor}.strelka2.indel.Somatic.vcf

	# compress temp folder
	tar -zcvf ${fileTumor}.strelka2.temp.tar.gz --remove-files ${fileTumor}.strelka2.temp
	check_file ${fileTumor}.strelka2.temp.tar.gz

	# mv files to own folder
	mkdir ${BaseName}_strelka2
	mv -f *.strelka2.* ./${BaseName}_strelka2/
	echo "*"
	echo "* END STRELKA"
	echo "*"
fi

######################################################################################################################################################################################################
### STEP 5: run bcftools norm to correct variant calling. Results will go to the {fileTumor}.bcftools.norm directory and will be re-annotated
######################################################################################################################################################################################################

if true; then
	date
	echo "*"
	echo "* BEGIN bcftools norm"
	echo "*"

	# for Mutect:
	date
	echo "** FOR MUTECT:"

	dir=${nameTumor}_mutect/

	/mnt/bioinfnas/prostate/daguilar/soft/strelka-2.9.2.centos6_x86_64/libexec/bgzip -c ${dir}${fileTumor}.Mutect.passed.vcf > ${dir}${fileTumor}.Mutect.passed.vcf.gz
	check_file ${dir}${fileTumor}.Mutect.passed.vcf.gz

	date
	echo "** Indexing vcf file"
	/mnt/bioinfnas/prostate/daguilar/soft/strelka-2.9.2.centos6_x86_64/libexec/tabix -p vcf ${dir}${fileTumor}.Mutect.passed.vcf.gz
	check_file ${dir}${fileTumor}.Mutect.passed.vcf.gz.tbi

	date
	echo "** Running bcftools"
	/mnt/bioinfnas/prostate/ldelgado/softwares/bcftools-1.17/bin/bcftools norm -f ${reference} -m -both -c w -o ${dir}${fileTumor}.Mutect.passed.vcf.norm -O v ${dir}${fileTumor}.Mutect.passed.vcf.gz
	check_file ${dir}${fileTumor}.Mutect.passed.vcf.norm

	
	# for varscan:
	date
	echo "** FOR VARSCAN:"
	
	dir=${nameTumor}_varscan

	for file in ${dir}/*.hc.vcf; do
		basename=${file%*.hc.vcf}
		/mnt/bioinfnas/prostate/daguilar/soft/strelka-2.9.2.centos6_x86_64/libexec/bgzip -c ${file} > ${basename}.hc.vcf.gz
	        check_file ${basename}.hc.vcf.gz
		
		date
		echo "** Indexing vcf file"
		/mnt/bioinfnas/prostate/daguilar/soft/strelka-2.9.2.centos6_x86_64/libexec/tabix -p vcf ${basename}.hc.vcf.gz
		check_file ${basename}.hc.vcf.gz.tbi
			
		date
		echo "** Running bcftools"
		/mnt/bioinfnas/prostate/ldelgado/softwares/bcftools-1.17/bin/bcftools norm -f ${reference} -m -both -c w -o ${basename}.hc.vcf.norm -O v ${basename}.hc.vcf.gz
		check_file ${basename}.hc.vcf.norm
	done


	# for Strelka2:
	date
	echo "** FOR STRELKA:"
	
	dir=${nameTumor}_strelka2/
	
	for file in ${dir}/*.vcf; do
                basename=${file%*.vcf}
		/mnt/bioinfnas/prostate/daguilar/soft/strelka-2.9.2.centos6_x86_64/libexec/bgzip -c ${basename}.vcf > ${basename}.vcf.gz
		check_file ${basename}.vcf.gz
				
		date
		echo "	** Indexing vcf file"
		/mnt/bioinfnas/prostate/daguilar/soft/strelka-2.9.2.centos6_x86_64/libexec/tabix -p vcf ${basename}.vcf.gz
		check_file ${basename}.vcf.gz.tbi

		date
		echo "  ** Running bcftools"
		/mnt/bioinfnas/prostate/ldelgado/softwares/bcftools-1.17/bin/bcftools norm -f ${reference} -m -both -c w -o ${basename}.vcf.norm -O v ${basename}.vcf.gz
		check_file ${basename}.vcf.norm
	done
fi

######################################################################################################################################################################################################
### STEP 6: copy number using VARSCAN (http://varscan.sourceforge.net/copy-number-calling.html)
######################################################################################################################################################################################################

if true; then
	date
	echo "*"
	echo "* BEGIN VARSCAN - COPY NUMBER"
	echo "*"

	# Run pileup
	# Varscan requires a pileup file from each of the bam files.
	echo "** Running pileup for both normal and tumoral bam files"
	samtools mpileup -C 50 -B -q 1 -Q 15 --positions ${intervals} -f ${reference} ${fileNormal} ${fileTumor} > ${fileNormal}.${fileTumor}.intervals.pileup
	check_file ${fileNormal}.${fileTumor}.intervals.pileup

	date
	echo "** Running VarScan copynumber"
	java -jar /mnt/bioinfnas/prostate/daguilar/soft/VarScan.v2.4.3.jar copynumber ${fileNormal}.${fileTumor}.intervals.pileup  --mpileup 1 
 	check_file output.copynumber

	mv output.copynumber ${nameTumor}.copynumber

	date
	echo "** Running VarScan copyCaller"
	java -jar /mnt/bioinfnas/prostate/daguilar/soft/VarScan.v2.4.3.jar copyCaller ${nameTumor}.copynumber --output-file ${nameTumor}.copynumber.called
	check_file ${nameTumor}.copynumber.called

	mv -f *.copynumber ./${BaseName}_varscan/
	mv -f *.copynumber.called ./${BaseName}_varscan/
	mv -f *.copynumber.called.gc ./${BaseName}_varscan/

	date
	echo "* END VARSCAN - COPY NUMBER"
fi

######################################################################################################################################################################################################
# mv source files to own folder
######################################################################################################################################################################################################

mkdir ${BaseName}_bams

mv -f ${BaseName}*.bam ./${BaseName}_bams/.
mv -f ${BaseName}*.bai ./${BaseName}_bams/.
mv -f ${BaseName}*.pileup ./${BaseName}_bams/.
#the files below were created by the script "WES.mapping.sh" but, since they sould be in the folder *_bams, here we move them
mv -f ${BaseName}*.bam.ValidateSamFile.txt ./${BaseName}_bams/.
mv -f ${BaseName}*.recalibration_data.table ./${BaseName}_bams/.
mv -f ${BaseName}*.sam.err ./${BaseName}_bams/.

echo "*** Finished successfully ***"
date
