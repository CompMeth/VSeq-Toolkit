outDir=$1
file1=$2
file2=$3
combinedRef=$4
stringencyCont=$5
UMthresholdCont=$6
aligner=$7
samtools=$8

cd $outDir
echo "Contaminant Analysis started" >> log
date >> log
echo " "
echo "Alignment in process..."
echo " "
echo "$aligner mem -H -M -t 8 $combinedRef $outDir/*pair1.fastq $outDir/*pair2.fastq > $outDir/alignment.sam"
$aligner mem -H -M -t 8 $combinedRef $outDir/*pair1.fastq $outDir/*pair2.fastq > $outDir/alignment.sam
$samtools view -bT $combinedRef $outDir/alignment.sam > $outDir/alignment.bam
$samtools sort $outDir/alignment.bam -o $outDir/alignment.sorted.bam
$samtools rmdup $outDir/alignment.sorted.bam  $outDir/alignment.nodup.bam
$samtools view $outDir/alignment.nodup.bam > $outDir/alignment.nodup.sam
awk '($2<'256')'  $outDir/alignment.nodup.sam > $outDir/alignment.nodup.remsec.sam
awk '($3!="*")' $outDir/alignment.nodup.remsec.sam > $outDir/alignment.nodup.remsec.alignedOnly.sam
echo " "
echo "Selection of reads in process..."
echo " "
if [[ $stringencyCont == *high* ]]; then
#	 echo "high"
	 awk '($2=="83" || $2=="163" || $2=="99" || $2=="147" )' $outDir/alignment.nodup.remsec.alignedOnly.sam | awk '($7=="=")' > $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam
fi
if [[ $stringencyCont == *medium* ]]; then
#	 echo "medium"
         awk '($2=="83" || $2=="163" || $2=="99" || $2=="147" || $2=="81" || $2=="161" || $2=="97" || $2=="145")' $outDir/alignment.nodup.remsec.alignedOnly.sam | awk '($7=="=")' > $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam
fi
if [[ $stringencyCont == *low* ]]; then
#	 echo "low"
         awk '($2=="83" || $2=="163" || $2=="99" || $2=="147" || $2=="81" || $2=="161" || $2=="97" || $2=="145" || $2=="65"  || $2=="129" || $2=="113" || $2=="177")' $outDir/alignment.nodup.remsec.alignedOnly.sam > $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam
fi

cut -f1 $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam | sort | uniq > $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam1
grep -Fwf $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam1 $outDir/alignment.sam > $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam2

if [[ $stringencyCont == *high* ]]; then
	awk '($2=="83" || $2=="163" || $2=="99" || $2=="147" )' $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam2 > $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam3
fi
if [[ $stringencyCont == *medium* ]]; then
         awk '($2=="83" || $2=="163" || $2=="99" || $2=="147" || $2=="81" || $2=="161" || $2=="97" || $2=="145")' $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam2 > $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam3
fi
if [[ $stringencyCont == *low* ]]; then
         awk '($2=="83" || $2=="163" || $2=="99" || $2=="147" || $2=="81" || $2=="161" || $2=="97" || $2=="145" || $2=="65"  || $2=="129" || $2=="113" || $2=="177")'  $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam2 > $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam3
fi
################

echo " "
echo "Distribution statistics generation in process..."
echo " "
awk '($2=="83" || $2=="99" ||  $2=="81" || $2=="97" || $2=="65" || $2=="113")' $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam3 | awk '{if (($3 ~ /_random/)  || ($3 ~ /chrUn/) || ($3 ~ /_alt/) ) {$3="Other-Chr"; print} else if  (($3 !~ /_random/)  || ($3 !~ /chrUn/) || ($3 !~ /_alt/) ) {print $0;}}' | sed 's/ /\t/g' | cut -f3 | sort | uniq -c > $outDir/temp1 

grep 'chr' $outDir/temp1 > $outDir/temp2
grep 'Chr' $outDir/temp1 >> $outDir/temp2
sed -i 's/^ *//g' $outDir/temp2


grep -v 'chr' $outDir/temp1 > $outDir/temp3
grep -v 'Chr' $outDir/temp3 > $outDir/temp4
sed -i 's/^ *//g' $outDir/temp4

awk '{print $2 "," $1}' $outDir/temp2 > $outDir/temp2a
echo  "chromsome""," "NumberOfReadPairs" > $outDir/header1
cat $outDir/header1 $outDir/temp2a > $outDir/ContAnalysis_NumberOfReadPairsPerChromsomes.csv
awk -F ',' '{sum+=$2}END{print "Chromsomes" "," sum}' $outDir/ContAnalysis_NumberOfReadPairsPerChromsomes.csv > f1

awk '{print $2 "," $1}' $outDir/temp4 > $outDir/temp4a
echo  "ContaminatOrReferenceSequence""," "NumberOfReadPairs" > $outDir/header2
cat $outDir/header2 $outDir/temp4a  $outDir/f1 > $outDir/ContAnalysis_NumberOfReadPairsPerContaminatORReferenceSequence.csv
 

sort -k1,1 -k2,2n  $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam3 | cut -f1,2,3,4,6,9,10,12,13,14 | sed 'N;s/\n/\t/g' | cut -f1,2,3,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20 > $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam4
awk '{if ($2=="99" || $2=="163" || $2=="97" || $2=="161" || $2=="65" || $2=="129" ) print $1 "\t" "+" "\t"  $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $19 ; else print $1 "\t" "-" "\t"  $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $19}' $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam4 | awk '{if ($11=="99" || $11=="163" || $2=="97" || $2=="161" || $2=="65" || $2=="129") print $1 "\t" $2 "\t"  $3 "\t" $4-1 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" "+" "\t" $12 "\t" $13-1 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $19 ; else print $1 "\t" $2 "\t"  $3 "\t" $4-1 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" "-" "\t" $12 "\t" $13-1 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $19}' | sed  's/AS:i://g' | sed  's/XS:i://g'  |  awk '{print $0 "\t" $10/$9 "\t" $19/$18}' | awk -v UMthresholdCont=$UMthresholdCont '{if ($20>=UMthresholdCont && $21>=UMthresholdCont) print $0 "\t" "NonUniqueAlignedReadPair"; else print $0 "\t" "UniqueAlignedReadPair"}' | sed 's/NM:i://g' > $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam5

echo "Read-ID,FirstReadStrand,FirstReadMappedToReference,FirstReadMappingPosition,FirstReadCIGAR-matchPartInfo,FirstReadInsertSizeOfReadPair,FirstReadSequence,FirstReadNumberOfMismatchBases,FirstReadPrimaryAlignmentScore, FirstReadSecondaryAlignmentScore,SecondReadStrand,SecondReadMappedToReference,SecondReadMappingPosition,SecondReadCIGAR-matchPartInfo,SecondReadInsertSizeOfReadPair,SecondReadSequence,SecondReadNumberOfMismatchBases,SecondReadPrimaryAlignmentScore, SecondReadSecondaryAlignmentScore,FirstReadRatioOfPrimarySecondaryRead,SecondReadRatioOfPrimarySecondaryRead,Reliability"  | cut -d ',' -f1-5,7,11-14,16,22 > $outDir/ContAnalysis_DetailsOfReadPairsPerChromsomes.csv
awk '$3  ~ /chr/' $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam5 | sed 's/\t/,/g ' | cut -d ',' -f1-5,7,11-14,16,22 >> $outDir/ContAnalysis_DetailsOfReadPairsPerChromsomes.csv
echo "Read-ID,FirstReadStrand,FirstReadMappedToReference,FirstReadMappingPosition,FirstReadCIGAR-matchPartInfo,FirstReadInsertSizeOfReadPair,FirstReadSequence,FirstReadNumberOfMismatchBases,FirstReadPrimaryAlignmentScore, FirstReadSecondaryAlignmentScore,SecondReadStrand,SecondReadMappedToReference,SecondReadMappingPosition,SecondReadCIGAR-matchPartInfo,SecondReadInsertSizeOfReadPair,SecondReadSequence,SecondReadNumberOfMismatchBases,SecondReadPrimaryAlignmentScore, SecondReadSecondaryAlignmentScore,FirstReadRatioOfPrimarySecondaryRead,SecondReadRatioOfPrimarySecondaryRead,Reliability" | cut -d ',' -f1-5,7,11-14,16,22  > $outDir/ContAnalysis_DetailsOfReadPairsPerContaminatSequence.csv
awk '$3 !~ /chr/' $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam5 | awk '$3 !~ /Chr/' |  sed 's/\t/,/g ' | grep -v 'vector' | cut -d ',' -f1-5,7,11-14,16,22 >>  $outDir/ContAnalysis_DetailsOfReadPairsPerContaminatSequence.csv


echo "Read-ID,FirstReadStrand,FirstReadMappedToReference,FirstReadMappingPosition,FirstReadCIGAR-matchPartInfo,FirstReadInsertSizeOfReadPair,FirstReadSequence,FirstReadNumberOfMismatchBases,FirstReadPrimaryAlignmentScore, FirstReadSecondaryAlignmentScore,SecondReadStrand,SecondReadMappedToReference,SecondReadMappingPosition,SecondReadCIGAR-matchPartInfo,SecondReadInsertSizeOfReadPair,SecondReadSequence,SecondReadNumberOfMismatchBases,SecondReadPrimaryAlignmentScore, SecondReadSecondaryAlignmentScore,FirstReadRatioOfPrimarySecondaryRead,SecondReadRatioOfPrimarySecondaryRead,Reliability" | cut -d ',' -f1-5,7,11-14,16,22 >  $outDir/ContAnalysis_DetailsOfReadPairsPerVectorReferenceSequence.csv
awk '$3 !~ /chr/' $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam5 | awk '$3 !~ /Chr/' |  sed 's/\t/,/g ' | grep 'vector' | cut -d ',' -f1-5,7,11-14,16,22 >>  $outDir/ContAnalysis_DetailsOfReadPairsPerVectorReferenceSequence.csv


echo " "
#echo "Combined final BAM file generation in process..."
#echo " "
#$samtools view -bS -t  $combinedRef.fai $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam3 > $outDir/ContAnalysis_CombinedCorrectlyAlignedReadPairs.bam
#echo " "
echo "FASTQ files without contaminants sequences are being generated..."
echo " "

awk '($3=="*")' $outDir/alignment.sam | cut -f1 | sort | uniq -d  > $outDir/notAligned
echo " "
echo "Approximate fragment size statistics generation in process..."
echo " "
#for chr
awk '($2=="83" || $2=="99" || $2=="97" || $2=="81" || $2=="65" || $2=="113" )' $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam3 | grep 'chr' | cut -f9  | sort -k1,1nr |  sed 's/-//g'  | awk '($1 <= 50)' | wc -l > $outDir/size
awk '($2=="83" || $2=="99" || $2=="97" || $2=="81" || $2=="65" || $2=="113" )' $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam3 | grep 'chr' | cut -f9  | sort -k1,1nr |  sed 's/-//g'  | awk '($1 > 50 && $1 <= 100)' | wc -l >> $outDir/size
awk '($2=="83" || $2=="99" || $2=="97" || $2=="81" || $2=="65" || $2=="113" )' $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam3 | grep 'chr' | cut -f9  | sort -k1,1nr |  sed 's/-//g'  | awk '($1 > 100 && $1 <= 200)' | wc -l >> $outDir/size
awk '($2=="83" || $2=="99" || $2=="97" || $2=="81" || $2=="65" || $2=="113" )' $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam3 | grep 'chr' | cut -f9  | sort -k1,1nr |  sed 's/-//g'  | awk '($1 > 200 && $1 <= 500)' | wc -l >> $outDir/size
awk '($2=="83" || $2=="99" || $2=="97" || $2=="81" || $2=="65" || $2=="113" )' $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam3 | grep 'chr' | cut -f9  | sort -k1,1nr |  sed 's/-//g'  |awk '($1 > 500 && $1 <= 1000)' | wc -l >> $outDir/size
awk '($2=="83" || $2=="99" || $2=="97" || $2=="81" || $2=="65" || $2=="113" )' $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam3 | grep 'chr' | cut -f9  | sort -k1,1nr |  sed 's/-//g'  | awk '($1 > 1000)' | wc -l >> $outDir/size


echo "size<=50", "size>50.<=100", "size>100.<=200", "size>200.<=500", "size>500.<=1000", "size>1000", > $outDir/header3
sed '$!{:a;N;s/\n/,/;ta}'  $outDir/size > $outDir/size1 
cat $outDir/header3 $outDir/size1 > $outDir/ContAnalysis_FragmentSizeDistribution.csv
rm $outDir/size $outDir/size1
#for all cont together
awk '($2=="83" || $2=="99" || $2=="97" || $2=="81" || $2=="65" || $2=="113" )' $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam3 | grep -v 'chr' | grep -v 'vector' | cut -f9  | sort -k1,1nr |  sed 's/-//g'  | awk '($1 <= 50)' | wc -l > $outDir/size
awk '($2=="83" || $2=="99" || $2=="97" || $2=="81" || $2=="65" || $2=="113" )' $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam3 | grep -v 'chr' | grep -v 'vector' | cut -f9  | sort -k1,1nr |  sed 's/-//g'  | awk '($1 > 50 && $1 <= 100)' | wc -l >> $outDir/size
awk '($2=="83" || $2=="99" || $2=="97" || $2=="81" || $2=="65" || $2=="113" )' $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam3 | grep -v 'chr' | grep -v 'vector' | cut -f9  | sort -k1,1nr |  sed 's/-//g'  | awk '($1 > 100 && $1 <= 200)' | wc -l >> $outDir/size
awk '($2=="83" || $2=="99" || $2=="97" || $2=="81" || $2=="65" || $2=="113" )' $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam3 | grep -v 'chr' | grep -v 'vector' | cut -f9  | sort -k1,1nr |  sed 's/-//g'  | awk '($1 > 200 && $1 <= 500)' | wc -l >> $outDir/size
awk '($2=="83" || $2=="99" || $2=="97" || $2=="81" || $2=="65" || $2=="113" )' $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam3 | grep -v 'chr' | grep -v 'vector' | cut -f9  | sort -k1,1nr |  sed 's/-//g'  |awk '($1 > 500 && $1 <= 1000)' | wc -l >> $outDir/size
awk '($2=="83" || $2=="99" || $2=="97" || $2=="81" || $2=="65" || $2=="113" )' $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam3 | grep -v 'chr' | grep -v 'vector' | cut -f9  | sort -k1,1nr |  sed 's/-//g'  | awk '($1 > 1000)' | wc -l >> $outDir/size

echo "size<=50", "size>50.<=100", "size>100.<=200", "size>200.<=500", "size>500.<=1000", "size>1000", > $outDir/header3
sed '$!{:a;N;s/\n/,/;ta}'  $outDir/size > $outDir/size1
less $outDir/size1 >> $outDir/ContAnalysis_FragmentSizeDistribution.csv
rm $outDir/size $outDir/size1

#for vector only
awk '($2=="83" || $2=="99" || $2=="97" || $2=="81" || $2=="65" || $2=="113" )' $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam3 | grep 'vector' | cut -f9  | sort -k1,1nr |  sed 's/-//g'  | awk '($1 <= 50)' | wc -l > $outDir/size
awk '($2=="83" || $2=="99" || $2=="97" || $2=="81" || $2=="65" || $2=="113" )' $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam3 | grep 'vector' | cut -f9  | sort -k1,1nr |  sed 's/-//g'  | awk '($1 > 50 && $1 <= 100)' | wc -l >> $outDir/size
awk '($2=="83" || $2=="99" || $2=="97" || $2=="81" || $2=="65" || $2=="113" )' $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam3 | grep 'vector' | cut -f9  | sort -k1,1nr |  sed 's/-//g'  | awk '($1 > 100 && $1 <= 200)' | wc -l >> $outDir/size
awk '($2=="83" || $2=="99" || $2=="97" || $2=="81" || $2=="65" || $2=="113" )' $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam3 | grep 'vector' | cut -f9  | sort -k1,1nr |  sed 's/-//g'  | awk '($1 > 200 && $1 <= 500)' | wc -l >> $outDir/size
awk '($2=="83" || $2=="99" || $2=="97" || $2=="81" || $2=="65" || $2=="113" )' $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam3 | grep 'vector' | cut -f9  | sort -k1,1nr |  sed 's/-//g'  |awk '($1 > 500 && $1 <= 1000)' | wc -l >> $outDir/size
awk '($2=="83" || $2=="99" || $2=="97" || $2=="81" || $2=="65" || $2=="113" )' $outDir/alignment.nodup.remsec.alignedOnly.CorrectlyPaired.sam3 | grep 'vector' | cut -f9  | sort -k1,1nr |  sed 's/-//g'  | awk '($1 > 1000)' | wc -l >> $outDir/size

echo "size<=50", "size>50.<=100", "size>100.<=200", "size>200.<=500", "size>500.<=1000", "size>1000", > $outDir/header3
sed '$!{:a;N;s/\n/,/;ta}'  $outDir/size > $outDir/size1
less $outDir/size1 >> $outDir/ContAnalysis_FragmentSizeDistribution.csv

echo "Reference,chr,contaminants_All,vector" > $outDir/head1
sed -i 's/,/\n/g' $outDir/head1
paste $outDir/head1 $outDir/ContAnalysis_FragmentSizeDistribution.csv | sed 's/\t/,/g' >> $outDir/ContAnalysis_FragmentSizeDistribution.csv.temp
mv $outDir/ContAnalysis_FragmentSizeDistribution.csv.temp $outDir/ContAnalysis_FragmentSizeDistribution.csv

#FASTQ fiel without conmtaminants seqeunces for further analysis
awk -F ','  '(($3 ~ /chr/ || $3 ~ /vector/) && ($8 ~ /chr/ || $8 ~ /vector/))'  $outDir/ContAnalysis_DetailsOfReadPairsPerChromsomes.csv  ContAnalysis_DetailsOfReadPairsPerVectorReferenceSequence.csv | cut -d ',' -f1 | sort | uniq | cat - notAligned > mappedToVecORChrORnotMapped
grep --no-group-separator -A3  -Fwf $outDir/mappedToVecORChrORnotMapped  $outDir/*pair1.fastq > $outDir/ContAnalysis_WithoutContaminantsReadPairs.R1.fastq
grep --no-group-separator -A3  -Fwf $outDir/mappedToVecORChrORnotMapped  $outDir/*pair2.fastq > $outDir/ContAnalysis_WithoutContaminantsReadPairs.R2.fastq
gzip $outDir/ContAnalysis_WithoutContaminantsReadPairs.R1.fastq $outDir/ContAnalysis_WithoutContaminantsReadPairs.R2.fastq
#mv $outDir/alignment.sam $outDir/alignment_ContAnalysis.sam
rm -f alignment.* header* f1 file* not* size* temp* head1 mappedToVecORChrORnotMapped *pair1.fastq *pair2.fastq
echo " "
echo "	  Contaminant Analysis Finished"
echo "Contaminant Analysis Finished" >> log
date >> log
echo "..........................................."
