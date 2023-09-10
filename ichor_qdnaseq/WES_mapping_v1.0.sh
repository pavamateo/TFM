#!/bin/bash
#SBATCH -p long
#SBATCH --job-name="WES_map"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=out.%x.%j
#SBATCH --error=err.%x.%j

### check the cpus to use. Using 20 cpus in the vhio cluster reduces the number of jobs you could run simultaneously in the server. Some samples with this features take more than 12hrs.  

### USER: this script has to be in the same folder as fastq files. Then, it all depends on how many fastq files are, since all the pipeline runs by using one pair of fastqR1 - fastqR2, and note that each step is going to be consecutively (not parallel).      
### USER in VHIO SERVER: sbatch WES_mapping_v1.sh 


RANDOM=$(date +%s%N | cut -b10-19)
RANDOMSTAMP='temp.'$RANDOM$RANDOM$RANDOM$RANDOM

### FOR ADAPT A PART OF THE CODE, WE CAN COMMENT A BLOCK OF CODE:
#: <<'END'
#.....CODE TO BE COMMENTED....
#END

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


### GENOME RELATED FILES

reference='/mnt/bioinfnas/prostate/daguilar/references/GRCh38_exomes/Homo_sapiens_assembly38.fasta'
knownsites='/mnt/bioinfnas/prostate/daguilar/references/GRCh38_exomes/Homo_sapiens_assembly38.dbsnp138.vcf'
kapa='/mnt/bioinfnas/prostate/daguilar/references/GRCh38_exomes/KAPA_HyperExome_primary_targets.bed'


check_file ${reference}
check_file ${knownsites}
check_file ${kapa}


### TRIMMING, FASTQC METRICS, MAPPING, QUALITY METRICS, MARKDUP AND RECALIBRATION 
for file in *_1.fastq.gz
do
	date

	basename=${file%*_1.fastq.gz}
	path=`readlink -f "${file}"`
        dirname=`dirname $path`

	f1=${basename}_1.fastq.gz
	f2=${basename}_2.fastq.gz

	if true; then

		date
                echo "** FASTQC"

		outdirqc=${dirname}/${basename}_fastqc
		mkdir ${outdirqc}

		/mnt/bioinfnas/prostate/daguilar/soft/FastQC/fastqc -o ${outdirqc} --threads 20 --nogroup ${f1} ${f2} --extract

		date
		echo "** Trimming"

		check_file "/mnt/bioinfnas/prostate/daguilar/references/illumina_adaptor_sequences.fa"

		java -jar /mnt/bioinfnas/prostate/daguilar/soft/Trimmomatic-0.39/trimmomatic-0.39.jar PE -phred33 -threads 20 ${f1} ${f2} ${basename}_paired_trimmed_R1.fastq ${basename}_unpaired_trimmed_R1.fastq ${basename}_paired_trimmed_R2.fastq ${basename}_unpaired_trimmed_R2.fastq ILLUMINACLIP:/mnt/bioinfnas/prostate/daguilar/references/illumina_adaptor_sequences.fa:2:30:10 MINLEN:110 &> ${basename}_trimmomatic.log

		pigz -v -f -p 25 ${basename}_paired_trimmed_R1.fastq
		pigz -v -f -p 25 ${basename}_paired_trimmed_R2.fastq
		pigz -v -f -p 25 ${basename}_unpaired_trimmed_R1.fastq
		pigz -v -f -p 25 ${basename}_unpaired_trimmed_R2.fastq

		R1=${basename}_paired_trimmed_R1.fastq.gz
		R2=${basename}_paired_trimmed_R2.fastq.gz

		check_file ${R1}
		check_file ${R2}
		
		date
		echo "** Mapping with BWA"

		/mnt/bioinfnas/prostate/daguilar/soft/bwa.kit/bwa mem -R "@RG\tID:1\tDS:HYPERCAP\tPL:ILLUMINA\tLB:SAMPLE\tSM:SAMPLE" -t 20 -M ${reference} ${R1} ${R2} | gzip > ${basename}.sam.gz 2> ${basename}.sam.err
		check_file ${basename}.sam.gz

		date
		echo "** SAM to BAM"

		/mnt/bioinfnas/prostate/daguilar/soft/gatk-4.2.5.0/gatk SamFormatConverter -I ${basename}.sam.gz -O ${basename}.nosort.bam --VERBOSITY ERROR --java-options " -Djava.io.tmpdir=/mnt/bioinfnas/prostate/ldelgado/GATK_tmp_files -Xmx32G"
		check_file ${basename}.nosort.bam

		date
                echo "*** gatk FixMateInformation"

                /mnt/bioinfnas/prostate/daguilar/soft/gatk-4.2.5.0/gatk FixMateInformation -I ${basename}.nosort.bam -O ${basename}.fixmate.bam --VERBOSITY ERROR --java-options " -Djava.io.tmpdir=/mnt/bioinfnas/prostate/ldelgado/GATK_tmp_files -Xmx32G"
                check_file ${basename}.fixmate.bam

		date
		echo "** Sorted BAM"
		
		/mnt/bioinfnas/prostate/daguilar/soft/gatk-4.2.5.0/gatk SortSam -I ${basename}.fixmate.bam -O ${basename}.bam -SO coordinate --VERBOSITY ERROR --java-options " -Djava.io.tmpdir=/mnt/bioinfnas/prostate/ldelgado/GATK_tmp_files -Xmx32G"
                check_file ${basename}.bam

		date
		echo "** Deleting SAM and preliminary BAM files to save disk space"

		rm -f ${basename}.sam.gz
		rm -f ${basename}.nosort.bam
		rm -f ${basename}.fixmate.bam

		date
		echo "** Verifying BAM files"

		java -Xmx32G -jar /mnt/bioinfnas/prostate/daguilar/soft/picard/picard.jar ValidateSamFile I=${basename}.bam OUTPUT=${basename}.bam.ValidateSamFile.txt MODE=VERBOSE MAX_RECORDS_IN_RAM=10000000

		if grep -q TRUNCATED_FILE ${basename}.bam.ValidateSamFile.txt; then
		    echo "$(date)	### ERROR: bam file ${basename}.bam is not valid. Check file ${basename}.bam.ValidateSamFile.txt";  exit
		fi
		
		date
		echo "** Marking duplicates"

		java  -Xmx32G  -jar /mnt/bioinfnas/prostate/daguilar/soft/picard/picard.jar MarkDuplicates VALIDATION_STRINGENCY=LENIENT INPUT=${basename}.bam OUTPUT=${basename}.markdups.bam SORTING_COLLECTION_SIZE_RATIO=0.10 MAX_RECORDS_IN_RAM=50000 METRICS_FILE=${basename}.picard_markduplicates_metrics.txt ASSUME_SORTED=true CREATE_INDEX=true
		check_file ${basename}.markdups.bam
		check_file ${basename}.markdups.bai
		rm -f ${basename}.bam

		date
		echo "** Recalibrating base quality scores step1"

		/mnt/bioinfnas/prostate/daguilar/soft/gatk-4.2.5.0/gatk BaseRecalibrator --input ${basename}.markdups.bam --known-sites ${knownsites} --reference ${reference} --output ${basename}.recalibration_data.table --java-options " -Djava.io.tmpdir=/mnt/bioinfnas/prostate/ldelgado/GATK_tmp_files -Xmx32G"
		check_file ${basename}.recalibration_data.table

		date
		echo "** Recalibrating base quality scores step2"
		/mnt/bioinfnas/prostate/daguilar/soft/gatk-4.2.5.0/gatk ApplyBQSR --bqsr-recal-file ${basename}.recalibration_data.table --input ${basename}.markdups.bam --output ${basename}.recal.bam --java-options " -Djava.io.tmpdir=/mnt/bioinfnas/prostate/ldelgado/GATK_tmp_files -Xmx32G"
		check_file ${basename}.recal.bam
		check_file ${basename}.recal.bai
		rm -f ${basename}.markdups.bam
		rm -f ${basename}.markdups.bai

		date
		echo "** CollectAlignmentSummaryMetrics"

		java -Xmx32G -jar /mnt/bioinfnas/prostate/daguilar/soft/picard/picard.jar CollectAlignmentSummaryMetrics METRIC_ACCUMULATION_LEVEL=ALL_READS I=${basename}.recal.bam O=${basename}.picard_alignment_metrics.txt R=${reference} VALIDATION_STRINGENCY=LENIENT
		check_file ${basename}.picard_alignment_metrics.txt

		date
		echo "** Running QUALIMAP"

		/mnt/bioinfnas/prostate/daguilar/soft/qualimap_v2.2.1/qualimap bamqc --java-mem-size=32G -bam ${basename}.recal.bam -outdir ${basename}_qualimap
		check_file ${basename}_qualimap/genome_results.txt

		date
		echo "** Count Reads"
		
		/mnt/bioinfnas/prostate/daguilar/soft/gatk-4.2.5.0/gatk CountReads -R ${reference} -I ${basename}.recal.bam -L ${kapa} -LE true > ${basename}.recal.bam.ontarget_reads.txt 
		check_file ${basename}.recal.bam.ontarget_reads.txt

		date
		echo "*** CollectInsertSizeMetrics"

	 	java -Xmx32G -jar /mnt/bioinfnas/prostate/daguilar/soft/picard/picard.jar CollectInsertSizeMetrics VALIDATION_STRINGENCY=LENIENT HISTOGRAM_FILE=${basename}.picard_insert_size_plot.pdf INPUT=${basename}.recal.bam OUTPUT=${basename}.picard_insert_size_metrics.txt
		check_file ${basename}.picard_insert_size_metrics.txt

		date
		echo "** CollectHsMetrics"
		/mnt/bioinfnas/prostate/daguilar/soft/gatk-4.2.5.0/gatk CollectHsMetrics --java-options -Xmx32G --BAIT_INTERVALS /mnt/bioinfnas/prostate/daguilar/references/GRCh38_exomes/KAPA_HyperExome_bait.interval_list --BAIT_SET_NAME KAPA_HyperExome --TARGET_INTERVALS /mnt/bioinfnas/prostate/daguilar/references/GRCh38_exomes/KAPA_HyperExome_target.interval_list --INPUT ${basename}.recal.bam --OUTPUT ${basename}.hs_metrics.txt --METRIC_ACCUMULATION_LEVEL ALL_READS --REFERENCE_SEQUENCE ${reference} -VALIDATION_STRINGENCY LENIENT --COVERAGE_CAP 100000 -PER_BASE_COVERAGE ${basename}.per_base_coverage.txt
		check_file ${basename}.hs_metrics.txt
		check_file ${basename}.per_base_coverage.txt
		
		date
		pigz -v -p 16 ${basename}.per_base_coverage.txt
		check_file ${basename}.per_base_coverage.txt.gz

		echo "** The BAM file from ${basename} is ready to be used"
		echo "*"

  	fi
	
	#save all the QC outputs in a single directory
	mkdir ${basename}_quality
	mv -f ${basename}.hs_metrics.txt ./${basename}_quality/.
	mv -f ${basename}.per_base_coverage.txt.gz ./${basename}_quality/.
	mv -f ${basename}.picard_* ./${basename}_quality/.
	mv -f ${basename}.recal.bam.ontarget_reads.txt ./${basename}_quality/.
	mv ${basename}_fastqc/ ./${basename}_quality/.
	mv ${basename}_qualimap/ ./${basename}_quality/.
	mv -f ${basename}_trimmomatic.log ./${basename}_quality/.

	#save the fastq files into a new folder for the processed fastqs
	mkdir ${basename}_trimmed
	mv -f ${basename}_trimmomatic.log ./${basename}_trimmed/.
	mv -f ${basename}_paired_trimmed* ./${basename}_trimmed/.
	mv -f ${basename}_unpaired_trimmed* ./${basename}_trimmed/.

	#save the renamed fastq files (the input files of this script) in a new folder
	mkdir ${basename}_fastqfiles
	mv -f ${f1} ./${basename}_fastqfiles/.
	mv -f ${f2} ./${basename}_fastqfiles/.

done

date
echo "*"
echo "*** Analysis finished successfully ***"
