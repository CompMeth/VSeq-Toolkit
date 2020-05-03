# README #
(VSeq-Toolkit 1.0)

VSeq-Toolkit is designed to analyze viral vector gene therapy sequencing data TES/WGS.

It has three modes Mode1;Contaminant Analysis, Mode2;Vector-Vector Fusion, Mode3;Vector-Host Fusion.

It is programmed in BASH and PERL.

### Setup of VSeq-Toolkit ###

## Installation procedure
Clone the repository from Github
```
git clone https://github.com/CompMeth/VSeq-Toolkit.git
cd VSeq-Toolkit

```
### Testing
To test if the installation of VSeq-Toolkit was successfull

```
First export the location of VSeq-toolkit
export VSeqToolkit=/complete-path/VSeq-Toolkit

Second run the test suite by executing following two commands on terminal
Execute the program within output directory
mkdir $VSeqToolkit/testDir/testResultsCheck
cd $VSeqToolkit/testDir/testResultsCheck
perl $VSeqToolkit/scripts/VSeq-TK.pl -c $VSeqToolkit/config.test.txt
OR
perl -I $VSeqToolkit/scripts/ $VSeqToolkit/scripts/VSeq-TK.pl -c $VSeqToolkit/config.test.txt

If process runs without error then installation is successfull 
The final result files in this $VSeqToolkit/testDir/testResults directory should be same as result files in $VSeqToolkit/testDir/testResultsCheck


```
## Dependencies

### Third-party tools
VSeq-Toolkit requires third-party tools that are already packaged within toolkit's 'thirdPartyTools' directory.
Tools (BWA-0.7.4, Skewer-0.1.117, Sametools-1.3.1, Bedtools-2.17.0)
The required Perl library is within the scripts directory.
```
$VSeqToolkit/thirdPartyTools
```

## Configuration file
The config.txt file can be used to analyze all three modes independently
OR Mode1, Mode2 and Mode3 
Only the relevant configuration file should be modified for particular analysis.
For testing GENE-IS installation user does not need to change any parameter in the configuration file.
The templates are in the gene-is path; i.e.
```
$VSeqToolkit/config.txt
```
VSeq-Toolkit only for research purpose

* Contact: saira.afzal@nct-heidelberg.de, saira.afzal@genewerk.de
