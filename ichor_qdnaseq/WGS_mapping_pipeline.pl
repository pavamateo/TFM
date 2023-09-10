

use strict;
use warnings;
use Getopt::Long;
use List::MoreUtils qw[any];
use Term::ANSIColor qw(:constants);
use File::Basename;
use Array::Utils qw(:all);
use List::MoreUtils qw[any uniq none distinct];
use List::Util qw[shuffle max min sum];
no if $] >= 5.017011, warnings => 'experimental::smartmatch';
no autovivification;

$Term::ANSIColor::AUTORESET = 1;

$| = 1;

use constant
	{
	READCOUNTER => "/home/daguilar/soft/hmmcopy_utils/bin/readCounter",
	REFERENCE => "/media/daguilar/SSD2TB/references/hg19",
	REFERENCE_FASTA => "/media/daguilar/SSD2TB/references/hg19/hg19_ref_genome.fasta",
	SAMTOOLS => "/usr/bin/samtools",
	BOWTIE2 => "/usr/bin/bowtie2",
	FASTQC => "/home/daguilar/soft/FastQC/fastqc",
	TRIMMOMATIC => "/home/daguilar/soft/Trimmomatic-0.39/trimmomatic-0.39.jar",
    TRIMGALORE => "/home/daguilar/soft/TrimGalore-0.6.7/trim_galore", 
	ILLUMINA_ADAPTOR_SEQUENCES => "/media/daguilar/SSD2TB/pipelines/WGS_CNV/illumina_adaptor_sequences.fa",
	THREADS => 12,
	STAMP => sprintf("%s%08s", "temp", int(rand(100000000))),
	TRUE => "1",
	FALSE => "0"
	};



my ($filename_pair1, $filename_pair2, $filename_output, $overwrite) = parse_arguments();

print "Started on: ".localtime()."\n\n";

print "Command-line arguments:\n";
print "\t--R1: $filename_pair1\n";
print "\t--R2: $filename_pair2\n";
print "\t--prefixes: $filename_output\n";
print "\t--overwrite: $overwrite\n";

print "\n";

print "Reading files:\n";

my @file_list1 = read_file($filename_pair1);
my @file_list2 = read_file($filename_pair2);
my @file_list3 = read_file($filename_output);

print "\n";

print "Checking data:\n";
sanity_check(\@file_list1, \@file_list2);
sanity_check(\@file_list1, \@file_list3);
print "\tLooks good.\n";
print "\n";

for my $i (0..$#file_list1)
	{
	print "Processing sample ".($i+1)." of ".scalar(@file_list1).":\n";

	my $file1 = $file_list1[$i];
	my $file2 = $file_list2[$i];
	my $output = $file_list3[$i];

	print "\tInput files:\n";
	print "\t\t$file1\n";
	print "\t\t$file2\n";
	print "\tOutput file prefix:\n";
	print "\t\t$output\n\n";

	die("Error! Output prefix '$output' looks like a folder, not a filename prefix!\n") if ($output =~ /\/$/);

	my ($name_output, $path_output) = fileparse($output);

	if ( ($overwrite == FALSE) and (-e "$output.500kb.wig") )
		{
		print "\t\tFile '$output.500kb.wig' already exists. ";
		print YELLOW, "Skipping.\n", RESET;
		}
	else
		{
        create_folder($path_output);

        # VERIFYPAIRED

        if (TRUE)
            {
	        system("/home/daguilar/soft/bbmap/reformat.sh in=$file1 in2=$file2 verifypaired &> $output.verifypaired.txt")
            }

        # FASTQC

        if (TRUE)
            {
	        system("mkdir ".$output."_fastqc") if (!-e $output."_fastqc");

	        system(FASTQC." $file1 -o ".${output}."_fastqc -t ".THREADS." --extract");
	        system(FASTQC." $file2 -o ".${output}."_fastqc -t ".THREADS." --extract");
            }



        print "\tTrimming sample ".($i+1)." of ".scalar(@file_list1).":\n";

        system_call("java -jar ".TRIMMOMATIC." PE -phred33 -threads ".THREADS." $file1 $file2 ".$output."_paired_trimmed_1.fastq ".$output."_unpaired_trimmed_1.fastq ".$output."_paired_trimmed_2.fastq ".$output."_unpaired_trimmed_2.fastq ILLUMINACLIP:".ILLUMINA_ADAPTOR_SEQUENCES.":2:30:10:8:FALSE -summary ".$output."_trimmomatic.log");


        system_call("pigz -fv -p ".THREADS." ".$output."_paired_trimmed_1.fastq", $output."_paired_trimmed_1.fastq.gz");
        system_call("pigz -fv -p ".THREADS." ".$output."_paired_trimmed_2.fastq", $output."_paired_trimmed_2.fastq.gz");
        system_call("pigz -fv -p ".THREADS." ".$output."_unpaired_trimmed_1.fastq", $output."_unpaired_trimmed_1.fastq.gz");
        system_call("pigz -fv -p ".THREADS." ".$output."_unpaired_trimmed_2.fastq", $output."_unpaired_trimmed_2.fastq.gz");

        $file1 = $output."_paired_trimmed_1.fastq.gz";
        $file2 = $output."_paired_trimmed_2.fastq.gz";
        
		print "\tMapping sample ".($i+1)." of ".scalar(@file_list1).":\n";

		my $string = BOWTIE2." -p ".THREADS." -x ".REFERENCE."/bowtie2_hg19 -1 $file1 -2 $file2 -S $output.bowtie2.sam  2> $output.bowtie2.log.txt";

		system("echo '$string' > $output.bowtie2.parameters.log.txt");

		system_call(BOWTIE2." -p ".THREADS." -x ".REFERENCE."/bowtie2_hg19 -1 $file1 -2 $file2 -S $output.bowtie2.sam  2> $output.bowtie2.log.txt", "$output.bowtie2.sam");

		system_call(SAMTOOLS." view -bS -@ ".THREADS." $output.bowtie2.sam > $output.bowtie2.bam", "$output.bowtie2.bam");

		system_call("rm $output.bowtie2.sam");

		system_call(SAMTOOLS." sort -@ ".THREADS." $output.bowtie2.bam -o $output.sorted.bowtie2.bam", "$output.sorted.bowtie2.bam");

		system_call("rm $output.bowtie2.bam");

		system_call(SAMTOOLS." index -@ ".THREADS." $output.sorted.bowtie2.bam");

		check_overall_alignment_rate("$output.bowtie2.log.txt");

		print "\tLOGs generation for sample ".($i+1)." of ".scalar(@file_list1).":\n";

		system_call(SAMTOOLS." sort -n -@ ".THREADS." -o ".$output."_namesort.bam $output.sorted.bowtie2.bam", $output."_namesort.bam");

		system_call(SAMTOOLS." fixmate -m -O bam -@ ".THREADS." ".$output."_namesort.bam ".$output."_fixmate.bam", $output."_fixmate.bam");

		system_call(SAMTOOLS." sort -@ ".THREADS." -o ".$output."_positionsort.bam ".$output."_fixmate.bam", $output."_positionsort.bam");

		system_call(SAMTOOLS." markdup -@ ".THREADS." ".$output."_positionsort.bam ".$output."_markdup.bam", $output."_markdup.bam");

		system_call(SAMTOOLS." rmdup -S ".$output."_markdup.bam ".$output.".rmdup.bam", $output.".rmdup.bam");

		system_call(SAMTOOLS." flagstat -@ ".THREADS." ".$output.".rmdup.bam > ".$output.".log");

		system_call(SAMTOOLS." index ".$output.".rmdup.bam");
		
		system_call("java -Xmx24g -jar /home/daguilar/soft/picard/picard.jar ValidateSamFile -I ".$output.".rmdup.bam -OUTPUT ".$output.".rmdup.bam.ValidateSamFile.txt -MODE SUMMARY -MAX_RECORDS_IN_RAM 10000000 -VERBOSITY ERROR 2>/dev/null", $output.".rmdup.bam.ValidateSamFile.txt");

		system_call("samtools flagstat --threads ".THREADS." $output.rmdup.bam > $output.rmdup.bam.flagstat.txt", "$output.rmdup.bam.flagstat.txt");

		system_call("/home/daguilar/soft/qualimap_v2.2.1/qualimap bamqc --java-mem-size=12G -bam ".$output.".recal.bam -outdir ".$output."_qualimap")

		system("rm -f  $output*namesort* ");
		system("rm -f  $output*fixmate* ");
		system("rm -f  $output*positionsort* ");

		print "\tCompressing temp .bam files to .cram files to save disk space:\n";

		system_call(SAMTOOLS." view -T ".REFERENCE_FASTA." -C -@ 20 -o ".$output."_markdup.cram ".$output."_markdup.bam", $output."_markdup.cram");

        system_call("rm -f ".$output."_markdup.bam");

		system_call(SAMTOOLS." view -T ".REFERENCE_FASTA." -C -@ 20 -o ".$output."_sorted.bowtie2.cram ".$output.".sorted.bowtie2.bam", $output."_sorted.bowtie2.cram");

        system_call("rm -f ".$output.".sorted.bowtie2.bam");



		print "\tGenerating WIG files:\n";

		system_call(READCOUNTER." --window 1000000 --quality 20 --chromosome 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY' $output.rmdup.bam > $output.1mb.wig", "$output.1mb.wig");
		system_call(READCOUNTER." --window 500000 --quality 20 --chromosome 'chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY' $output.rmdup.bam > $output.500kb.wig", "$output.500kb.wig");
		}
	}


print "Ended on: ".localtime()."\n\n";

####################################################################################################

sub check_file
	{
	my $file = $_[0];

	die("ERROR in file '$file'. EXITING.\n") if ( (-z $file) or (!-e $file) );
	}

sub parse_arguments
	{
	my ($filename_pair1, $filename_pair2, $filename_output, $overwrite);

	GetOptions("R1=s"=> \$filename_pair1,   
            "R2=s"  => \$filename_pair2,     
            "prefixes=s" => \$filename_output,
            "overwrite=s" => \$overwrite)
	or die "USAGE: $0 --R1 --R2 --prefixes [--overwrite]\n
\t--R1 = text file with the list of R1 files (with full path)
\t--R2 = text file with the list of R2 files (with full path)
\t--prefixes = text file with the list of prefixes for the output files (with full path)\n";

	if ( (!defined($filename_pair1)) or  (!defined($filename_pair2)) or  (!defined($filename_output)) )
		{
		die "USAGE: $0 --R1 --R2 --prefixes [--overwrite]\n
\t--R1 = text file with the list of R1 files (with full path)
\t--R2 = text file with the list of R2 files (with full path)
\t--prefixes = text file with the list of prefixes for the output files (with full path)\n";
		}

	if (!defined($overwrite))
		{$overwrite = FALSE;}
	elsif ($overwrite =~ /^true$/i)
		{$overwrite = TRUE;}
	elsif ($overwrite =~ /^false$/i)
		{$overwrite = FALSE;}
	else
		{
		die "USAGE: $0 --R1 --R2 --prefixes [--overwrite]\n
\t--R1 = text file with the list of R1 files (with full path)
\t--R2 = text file with the list of R2 files (with full path)
\t--prefixes = text file with the list of prefixes for the output files (with full path)\n";
		}

	return($filename_pair1, $filename_pair2, $filename_output, $overwrite);
	}


sub check_overall_alignment_rate
	{
	my $file = $_[0];

	my $rate = `grep 'overall alignment rate' $file`;

	($rate) = $rate =~ /(\S+)%/;

	if (!defined($rate))
		{
		print RED, "\n\t\tOverall alignment rate cannot be read!\n\n", RESET;
		}
	elsif ($rate > 80)
		{
		print GREEN, "\n\t\tOverall alignment rate: $rate%\n\n", RESET;
		}
	else
		{
		print RED, "\n\t\tOverall alignment rate: $rate%\n\n", RESET;
		}
	}

sub system_call
	{
	my ($string, $output_file) = @_;

	print "\t\tSystem call: ";
	print YELLOW, "$string\n", RESET;

	system($string);

	die("Error!\n") if ( (defined($output_file)) and (!-e $output_file) );
	}

sub create_folder
	{
	my $path_mapped = $_[0];

	if (!-e $path_mapped)
		{
		unless(mkdir $path_mapped)
			{die "Error! Unable to create folder '$path_mapped'\n";}
		}
	}


sub sanity_check
	{
	my ($file_list1_ref, $file_list2_ref) = @_;

	if (scalar(@$file_list1_ref) != scalar(@$file_list2_ref) )
		{die("Error! '$filename_pair1' and '$filename_pair2' have a different number of files! Are they really paired?\n");}

	foreach my $file (@$file_list1_ref)
		{
		die("Error! File '$file' has an incorrect name!\n") if ($file !~ /\w/);
		}

	foreach my $file (@$file_list2_ref)
		{
		die("Error! File '$file' has an incorrect name!\n") if ($file !~ /\w/);
		}
	}


sub read_file
    {
    my $file = $_[0];

	print "\tReading file '$file'... ";

	my @list;

    open(FHANDLE, $file) or die("\nError! File '$file' not found!\n");
	while(<FHANDLE>)
		{
		chomp;

		push(@list, $_);
		}
	close(FHANDLE);

	print "Done.\n";

	if (scalar(@list) == 0)
		{die("Error! '$file' seems to be empty!\n");}

	print "\tNumber of files in '$file': ".scalar(@list)."\n";
	return(@list);
    }
