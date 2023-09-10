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
	ICHORCNAR => "/media/daguilar/SSD2TB/pipelines/WGS_CNV/ichorCNA/runIchorCNA.R",
	GENECALLS => "/media/daguilar/SSD2TB/pipelines/WGS_CNV/ichorCNA/gene_calls_corrected_by_TF7.R",
    PANEL => "/media/daguilar/SSD2TB/references/panel.v3.txt",
    PANELREFFLAT => "/media/daguilar/SSD2TB/references/hg19/refFlat.panel.RData",
	STAMP => sprintf("%s%08s", "temp", int(rand(100000000))),
	TRUE => "1",
	FALSE => "0"
	};

# para ejecutar en el directorio que contiene los ficheros .wig
# este script ejecuta: gene_calls_corrected_by_TF.R

# perl /home/daguilar/pipelines/WGS_CNV/ichorCNA/ichorCNA_withPoN.pl

die("ERROR! No '".ICHORCNAR."'!\n") if (!-e ICHORCNAR);
die("ERROR! No '".GENECALLS."'!\n") if (!-e GENECALLS);
die("ERROR! No '".PANEL."'!\n") if (!-e PANEL);
die("ERROR! No '".PANELREFFLAT."'!\n") if (!-e PANELREFFLAT);

my $string = $ARGV[0];

die("ERROR! Must supply a file search string (e.g. 'PRO*FFPE*.500kb.wig') as argument (within quotes please!)\n") if (!defined($string));

print "Started on: ".localtime()."\n\n";

my @files = glob $string;

print "File list:\n";
foreach my $file (@files)
	{
	print "\t$file\n";
	}
print "\n";

my $counter = 1;

foreach my $file (@files)
	{
	my ($name, $path) = fileparse($file);

	my ($outdir_defaultPON) = $name =~ /(\S+).wig$/;

	$outdir_defaultPON = $path.$outdir_defaultPON."_default_PON_ichorCNA";

	my ($outdir_50salPON) = $name =~ /(\S+).wig$/;

	$outdir_50salPON = $path.$outdir_50salPON."_50salPON_ichorCNA";

	print "Sample $counter of ".scalar(@files)." (file search string: $string)\n";
	print "Input file:\n";
	print "\t$file\n";



    ###########################################
    # TRUE para el normal panel default de ichorCNA
    ###########################################

    if (TRUE)
        {
	    print "Output file:\n";
	    print "\t$outdir_defaultPON\n";

	    if (!defined($outdir_defaultPON))
		    {
		    die("Error! Cannot parse filename '$file'! Is it a .wig file?\n");
		    }

	    if (-e $outdir_defaultPON)
			    {
			    print "\t\tFolder '$outdir_defaultPON' already exists.\n";
			    print YELLOW, "\t\tSkipping.\n\n", RESET;
			    }
	    else
			    {
			    print "\tStarted on: ".localtime()."\n";
			    my $gcWig;
			    my $mapWig;
			    my $normalPanel;

			    if ($file =~ /\.500kb\./)
				    {
				    $gcWig = "/media/daguilar/SSD2TB/references/extdata/gc_hg19_500kb.wig";
				    $mapWig = "/media/daguilar/SSD2TB/references/extdata/map_hg19_500kb.wig";
				    $normalPanel = "/media/daguilar/SSD2TB/references/PONs/hg19/ichorCNA/HD_ULP_PoN_500kb_median_normAutosome_mapScoreFiltered_median.rds";
				    }
			    elsif ($file =~ /\.1mb\./)
				    {
				    $gcWig = "/media/daguilar/SSD2TB/references/extdata/gc_hg19_1000kb.wig";
				    $mapWig = "/media/daguilar/SSD2TB/references/extdata/map_hg19_1000kb.wig";
				    $normalPanel = "/media/daguilar/SSD2TB/references/PONs/hg19/ichorCNA/HD_ULP_PoN_1Mb_median_normAutosome_mapScoreFiltered_median.rds";
				    }
			    else
				    {die("ERROR! gcWig file undefined!\n");}

			    die("ERROR! No $gcWig!\n") if (!-e $gcWig);
			    die("ERROR! No $mapWig!\n") if (!-e $mapWig);
			    die("ERROR! No $normalPanel!\n") if (!-e $normalPanel);


			    create_folder($outdir_defaultPON);



			    my $string = "Rscript ".ICHORCNAR. " --id $name --WIG $file --ploidy 'c(2,3,4,6,8)' --normal 'c(0,0.5,0.6,0.7,0.8,0.9,0.95, 0.99, 0.995, 0.999)' --maxCN 5 --gcWig $gcWig --mapWig $mapWig --centromere /media/daguilar/SSD2TB/references/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt --normalPanel $normalPanel --includeHOMD True --chrTrain 'c(1:22)' --estimateNormal True --estimatePloidy True --estimateScPrevalence True --scStates 'c(1,3)' --txnE 0.9999 --txnStrength 10000 --normalizeMaleX True --outDir $outdir_defaultPON\n";
		    
			    system_call($string);       
                }

            if (FALSE)
                {
                $string = "Rscript ".GENECALLS." $outdir_50salPON/$file ".PANEL." ".PANELREFFLAT;

			    system_call($string);

                if (!-e "./$outdir_defaultPON/$file.gene.calls.panel.tsv")
                    {die("ERROR! No ./$outdir_defaultPON/$file.gene.calls.panel.tsv!\n")}

			    print "\n\tFinished on: ".localtime()."\n";
			    }
        }


    ###########################################
    # TRUE para el normal panel de 50 salivas
    ###########################################

    if (FALSE)
        {
	    print "Output file:\n";
    	print "\t$outdir_50salPON\n";

	    if (!defined($outdir_50salPON))
		    {
		    die("Error! Cannot parse filename '$file'! Is it a .wig file?\n");
		    }

	    if (-e $outdir_50salPON)
			    {
			    print "\t\tFolder '$outdir_50salPON' already exists.\n";
			    print YELLOW, "\t\tSkipping.\n\n", RESET;
			    }
	    else
			    {
			    print "\tStarted on: ".localtime()."\n";
			    my $gcWig;
			    my $mapWig;
			    my $normalPanel;

			    if ($file =~ /\.500kb\./)
				    {
				    $gcWig = "/media/daguilar/SSD2TB/references/extdata/gc_hg19_500kb.wig";
				    $mapWig = "/media/daguilar/SSD2TB/references/extdata/map_hg19_500kb.wig";
				    $normalPanel = "/media/daguilar/SSD2TB/references/PONs/hg19/ichorCNA/50sal_500kb_PON_median.rds";
                    
				    }
			    elsif ($file =~ /\.1mb\./)
				    {
				    $gcWig = "/media/daguilar/SSD2TB/references/extdata/gc_hg19_1000kb.wig";
				    $mapWig = "/media/daguilar/SSD2TB/references/extdata/map_hg19_1000kb.wig";
				    $normalPanel = "/media/daguilar/SSD2TB/references/PONs/hg19/ichorCNA/50sal_1000kb_PON_median.rds";
				    }
			    else
				    {die("ERROR! gcWig file undefined!\n");}

			    die("ERROR! No $gcWig!\n") if (!-e $gcWig);
			    die("ERROR! No $mapWig!\n") if (!-e $mapWig);
			    die("ERROR! No $normalPanel!\n") if (!-e $normalPanel);


			    create_folder($outdir_50salPON);

			    my $string = "Rscript ".ICHORCNAR. " --id $name --WIG $file --ploidy 'c(2,3,4,6,8)' --normal 'c(0,0.5,0.6,0.7,0.8,0.9,0.95, 0.99, 0.995, 0.999)' --maxCN 5 --gcWig $gcWig --mapWig $mapWig --centromere /media/daguilar/SSD2TB/references/extdata/GRCh37.p13_centromere_UCSC-gapTable.txt --normalPanel $normalPanel --includeHOMD True --chrTrain 'c(1:22)' --estimateNormal True --estimatePloidy True --estimateScPrevalence True --scStates 'c(1,3)' --txnE 0.9999 --txnStrength 10000 --normalizeMaleX True --outDir $outdir_50salPON\n";
		    
			    system_call($string);

                }

        if (FALSE)
                {
                $string = "Rscript ".GENECALLS." $outdir_50salPON/$file ".PANEL." ".PANELREFFLAT;

			    system_call($string);

                if (!-e "./$outdir_50salPON/$file.gene.calls.panel.tsv")
                    {die("ERROR! No ./$outdir_50salPON/$file.gene.calls.panel.tsv!\n")}
                }
        print "\n\tFinished on: ".localtime()."\n";
			      
        }
	$counter ++;
	}

print "Finished on: ".localtime()."\n\n";

####################################################################################################


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
			{die "Error! Unable to create folder '$path_mapped'\nDid you forget to create the './ichorCNA/' folder first?\n";}
		}
	}


sub parse_arguments
	{
	my ($filename_wig, $overwrite);


	GetOptions("wig=s"=> \$filename_wig,   
            "overwrite=s" => \$overwrite)
	or die "USAGE: $0 --wig [--overwrite]\n
\t--wig = path to .wig file
\n";


	if (!defined($filename_wig)) 
		{
		die "USAGE: $0 --wig [--overwrite]\n
\t--wig = path to .wig file
\n";
		}

	if (!defined($overwrite))
		{$overwrite = FALSE;}
	elsif ($overwrite =~ /^true$/i)
		{$overwrite = TRUE;}
	elsif ($overwrite =~ /^false$/i)
		{$overwrite = FALSE;}
	else
		{
		die "USAGE: $0 --wig [--overwrite]\n
\t--wig = path to .wig file
\n";
		}

	return($filename_wig, $overwrite);
	}



sub read_file
    {
    my $file = $_[0];

	print "Reading files in: $file... ";

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
