#!/usr/bin/perl -w
###################################################################################
###################################################################################
#										  #
#   VSeq-Toolkit- An automated toolkit for viral vector analyses                  #
#	1. Contmainants Distribution Analysis					  #
#	2. Vector-vector Fusion Analyis                         		  #
#	3. Vector-Host Fusion Analsyis	         				  #
#										  #
###################################################################################
###################################################################################
#										  #
#   VSeq-Toolkit 							          #
#   Version 1.0								          #
#   Author: Saira Afzal
#   Last update: 22-Dec-2019 							  #
#   VSeq-Tk.pl is the main interfacte program					  #
#   perl VSeq-Tk.pl -c config.txt						  #
#										  #
###################################################################################
###################################################################################
use Getopt::Std;
use Getopt::Long;
use Cwd;
use FindBin;
use lib "$FindBin::Bin";
use GeneisLib;

my @usage;
push @usage, "VSeq-TK: A toolkit for viral analysis.\n";
push @usage, "Version: 1.0 \n\n";
push @usage, "Usage: VSeq-TK.pl -c <configuration file> [options]\n\n";
push @usage, "Options:\n";
push @usage, "  -h, --help     Displays this information\n";
push @usage, "  -c, --configFile   Configuration file <required>\n";


my $help;
my $configFile;

GetOptions
(
 'h|help|?'    => \$help,
 'configFile=s'    => \$configFile,
);

if ($help) {
   print @usage;
   exit(0);
}
if (defined $configFile) {
   if (!-e $configFile){
        print "\nThe required configuration file $configFile does not exist!\n\n";
        print @usage;
        exit;
   }
}else{
    print "Please provide the required configuration file. Consult manual for details !!!\n\n";
    print @usage;
    exit;
}

#my $config = $configFile
my $config = new();	
$config->read($configFile);
my $outDir = $config->get_value("outDir");
my $bin = $config->get_value("bin");
my $trimmer = $config->get_value("trimmer");
my $aligner = $config->get_value("aligner");
my $samtools = $config->get_value("samtools");
my $mode = $config->get_value("mode"); 

my $file1 = $config->get_value("file1");
my $file2 = $config->get_value("file2");
my $qua = $config->get_value("qua");
my $lenPer = $config->get_value("lenPer");
my $adapter1 = $config->get_value("adapter1");
my $adapter2 = $config->get_value("adapter2");

system "echo  ..............................................";
print  "\n\n   QUALITY CONTROL \n\n";

system("echo bash  $bin/qualControl.sh $outDir $file1 $file2 $qua $lenPer $trimmer ");
system(" bash  $bin/qualControl.sh $outDir $file1 $file2 $qua $lenPer $adapter1 $adapter2 $trimmer ");


system "echo  ..............................................\n\n";
if ($config->has_value("contAna")){
	$contAna=$config->get_value("contAna");
	if ($contAna eq "true") {
		$combinedRef=$config->get_value("combinedRef");
		$stringencyCont=$config->get_value("stringencyCont");
		$UMthresholdCont=$config->get_value("UMthresholdCont");
				
		print "\n\n CONTAMINANTS ANALYSIS \n\n";
		system("echo bash  $bin/contAna.sh $outDir $outDir/file-pair1.fastq $outDir/file-pair2.fastq $combinedRef $stringencyCont $UMthresholdCont $aligner $samtools");
		system(" bash  $bin/contAna.sh $outDir $outDir/file-pair1.fastq $outDir/file-pair2.fastq $combinedRef $stringencyCont $UMthresholdCont $aligner $samtools");
	};
};

system "echo  ..............................................";

if ($config->has_value("vecVecFusion")){
        $vecVecFusion=$config->get_value("vecVecFusion");
        if ($vecVecFusion eq "true") {
                $vecRef=$config->get_value("vecRef");
                $stringencyVec=$config->get_value("stringencyVec");
                $UMthresholdVec=$config->get_value("UMthresholdVec");
                $minMapSpanVec=$config->get_value("minMapSpanVec");
		$distVecVec=$config->get_value("distVecVec");
		$opVecVec=$config->get_value("opVecVec");
		$idenVecVec=$config->get_value("idenVecVec");
		$contCheck=$config->get_value("contAna");

                print "\n\n VECTOR-VECTOR FUSION ANALYSIS \n\n";
                system("echo bash  $bin/vecVecFusion.sh $outDir $outDir/file-pair1.fastq $outDir/file-pair2.fastq $vecRef $stringencyVec $UMthresholdVec $minMapSpanVec $aligner $samtools $distVecVec $opVecVec $idenVecVec $bin $mode $contCheck");
                system("  bash  $bin/vecVecFusion.sh $outDir $outDir/file-pair1.fastq $outDir/file-pair2.fastq $vecRef $stringencyVec $UMthresholdVec $minMapSpanVec $aligner $samtools $distVecVec $opVecVec $idenVecVec $bin $mode $contCheck");
        };
};


system "echo  ..............................................";

if ($config->has_value("vecGenIS")){
        $vecGenIS=$config->get_value("vecGenIS");
        if ($vecGenIS eq "true") {
                $vecRef=$config->get_value("vecRef");
                $stringencyVecGen=$config->get_value("stringencyVecGen");
                $UMthresholdVecGen=$config->get_value("UMthresholdVecGen");
                $minMapSpanVecGen=$config->get_value("minMapSpanVecGen");
		$vecGenRef=$config->get_value("vecGenRef");
		$distVecGen=$config->get_value("distVecGen");
		$opVecGen=$config->get_value("opVecGen");
		$idenVecGen=$config->get_value("idenVecGen");
		$clusterRange=$config->get_value("clusterRange");
		$bedtools=$config->get_value("bedtools");
		$annoTable=$config->get_value("annoTable");
                $bin=$config->get_value("bin");
		$vecVecFusionCheck=$config->get_value("vecVecFusion");
		$contCheck=$config->get_value("contAna");		

                print "\n\n VECTOR-GENOME FUSION ANALYSIS \n\n";
                system("echo bash  $bin/vecGenIS.sh $outDir $outDir/file-pair1.fastq $outDir/file-pair2.fastq $vecRef $stringencyVecGen $UMthresholdVecGen $minMapSpanVecGen $aligner $samtools $vecGenRef $distVecGen $opVecGen $idenVecGen $clusterRange $bedtools $annoTable $bin $vecVecFusionCheck $mode $contCheck");
                system(" bash  $bin/vecGenIS.sh $outDir $outDir/file-pair1.fastq $outDir/file-pair2.fastq $vecRef $stringencyVecGen $UMthresholdVecGen $minMapSpanVecGen $aligner $samtools $vecGenRef $distVecGen $opVecGen $idenVecGen $clusterRange $bedtools $annoTable $bin $vecVecFusionCheck $mode $contCheck");
        };
};

system "echo  ..............................................";
