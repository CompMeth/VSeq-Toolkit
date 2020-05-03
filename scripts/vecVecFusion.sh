outDir=$1
file1=$2
file2=$3
vecRef=$4
stringencyVec=$5
UMthresholdVec=$6
minMapSpanVec=$7
aligner=$8
samtools=$9
distVecVec=${10}
opVecVec=${11}
idenVecVec=${12}
bin=${13}
mode=${14}
contCheck=${15}

cd $outDir
echo "Vector-Vector Fusion Analysis started" >> log
date >> log
echo " "
echo "Alignment with viral vector reference in process..."
echo " "

echo " "
if [[ $contCheck == *true* ]]; then
        $aligner mem -M -t 10 -T 15 $vecRef $outDir/ContAnalysis_WithoutContaminantsReadPairs.R1.fastq.gz  $outDir/ContAnalysis_WithoutContaminantsReadPairs.R2.fastq.gz > $outDir/alignmentV.sam
fi

if [[ $contCheck == *false* ]]; then
	$aligner mem -M -t 10 -T 15  $vecRef $outDir/*pair1.fastq $outDir/*pair2.fastq > $outDir/alignmentV.sam
fi
#Sensitive mode
if [[ $mode == *sensitive* ]]; then 
	echo ""
	echo "Sensitive mode selected"
	echo "Processing started..."
	mv alignmentV.sam alignmentSens.sam
	bash $bin/sensModule.sh $outDir $vecRef $aligner
	mv alignmentSens.sam.1 alignmentV.sam
	mv alignmentSens.3.sam.ext alignmentV.3.sam.ext
fi

mv $outDir/alignmentV.sam $outDir/alignmentV.3.sam
	
echo " "
echo "Selection of potential fusion reads in process..."
if [[ $stringencyVec == *high* ]]; then
        awk '($2=="99" || $2=="147" || $2=="83" || $2=="163")' $outDir/alignmentV.3.sam  | grep 'XP:Z:' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 | sed 's/XP:Z://g' > $outDir/alignmentV.4.sam
fi
if [[ $stringencyVec == *medium* ]]; then
        awk '($2=="99" || $2=="147" || $2=="83" || $2=="163" || $2=="81" || $2=="161" || $2=="97" || $2=="145")' $outDir/alignmentV.3.sam  | grep 'XP:Z:' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 | sed 's/XP:Z://g' > $outDir/alignmentV.4.sam
fi
if [[ $stringencyVec == *low* ]]; then
	awk '($2=="99" || $2=="147" || $2=="83" || $2=="163" || $2=="81" || $2=="161" || $2=="97" || $2=="145"  || $2=="65"  || $2=="129" || $2=="113" || $2=="177")' $outDir/alignmentV.3.sam  | grep 'XP:Z:' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 | sed 's/XP:Z://g' > $outDir/alignmentV.4.sam
fi
if [[ $stringencyVec == *null* ]]; then
        awk '($2=="99" || $2=="147" || $2=="83" || $2=="163" || $2=="81" || $2=="161" || $2=="97" || $2=="145"  || $2=="65"  || $2=="129" || $2=="113" || $2=="177" || $2=="131" || $2=="387" )' $outDir/alignmentV.3.sam  | grep 'XP:Z:' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 | sed 's/XP:Z://g' > $outDir/alignmentV.4.sam
fi

echo " "
echo "Post-processing of extracted reads..."
echo " "
cat $outDir/alignmentV.4.sam  |  cut -f15 | cut -d ';' -f1 | paste $outDir/alignmentV.4.sam - | cut -f1-14,16  | awk '{print $0";"}' > $outDir/alignmentV.4.sam.temp00
cut -f15  $outDir/alignmentV.4.sam.temp00|  awk '{print gsub (/;/, "")}' | paste $outDir/alignmentV.4.sam.temp00 - | awk '($16>=2)' | cut -f1 | grep -v -Fwf - $outDir/alignmentV.4.sam.temp00 > $outDir/alignmentV.4.sam.temp1
mv $outDir/alignmentV.4.sam.temp00 $outDir/alignmentV.4.sam

cut -f15 $outDir/alignmentV.4.sam  | sed 's/,/\t/g' > $outDir/alignmentV.5.sam

paste $outDir/alignmentV.4.sam $outDir/alignmentV.5.sam | awk '($3==$16)' > $outDir/alignmentV.6a.sam
######################

awk '{print $0, gsub(/S/,"\t",$6) }' alignmentV.6a.sam | awk '{print $0, gsub(/S/,"\t",$18) }' | sed 's/ /\t/g' > alignmentV.6a.sam.A

awk '($21==1)' alignmentV.6a.sam.A | awk '{ split($6,a , "S"); print $0 "\t" a[1];}' > alignmentV.6a.sam.B
awk '($23 ~ /M/ )' alignmentV.6a.sam.B | cut -f1-22 > alignmentV.6a.sam.B-MbeforeS
awk '($23 !~ /M/ )' alignmentV.6a.sam.B | cut -f1-22 > alignmentV.6a.sam.B-SbeforeM
awk '($21==2)' alignmentV.6a.sam.A  > alignmentV.6a.sam.C-2S

awk '{print $6}' alignmentV.6a.sam.B-MbeforeS | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*S//g' -e 's/[A-Z]//g'  | awk '{tot=0; for(i=1;i<=NF;i++) tot+=$i; print tot}' | awk '{print "0" "\t" $1}' | paste alignmentV.6a.sam.B-MbeforeS - > alignmentV.6a.sam.B-MbeforeS.1
awk '{print $6}' alignmentV.6a.sam.B-MbeforeS.1 | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*I//g' -e 's/[0-9]*M//g' -e 's/[A-Z]//g' | paste alignmentV.6a.sam.B-MbeforeS.1 - > alignmentV.6a.sam.B-MbeforeS.2

 awk '{print $6}' alignmentV.6a.sam.B-SbeforeM | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*I//g' -e 's/[0-9]*M//g' -e 's/[A-Z]//g' | paste alignmentV.6a.sam.B-SbeforeM - > alignmentV.6a.sam.B-SbeforeM.1
 awk '{print $6}' alignmentV.6a.sam.B-SbeforeM.1 | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*S//g' -e 's/[A-Z]//g'  | awk '{tot=0; for(i=1;i<=NF;i++) tot+=$i; print tot}' | awk '{print $1 "\t" "0"}' | paste alignmentV.6a.sam.B-SbeforeM.1 - > alignmentV.6a.sam.B-SbeforeM.2 

 awk '{print $6}' alignmentV.6a.sam.C-2S | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*I//g' -e 's/[0-9]*M//g' -e 's/[A-Z]//g' |  awk '{print $1 }' | paste alignmentV.6a.sam.C-2S - > alignmentV.6a.sam.C-2S.1
awk '{print $6}' alignmentV.6a.sam.C-2S.1 | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*S//g' -e 's/[A-Z]//g'  | awk '{tot=0; for(i=1;i<=NF;i++) tot+=$i; print tot}' | awk '{print $1 }' | paste alignmentV.6a.sam.C-2S.1 - > alignmentV.6a.sam.C-2S.2a
 awk '{print $6}'  alignmentV.6a.sam.C-2S.2a | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*I//g' -e 's/[0-9]*M//g' -e 's/[A-Z]//g' |  awk '{print $2 }' | paste alignmentV.6a.sam.C-2S.2a - > alignmentV.6a.sam.C-2S.2
cat alignmentV.6a.sam.B-MbeforeS.2  alignmentV.6a.sam.B-SbeforeM.2 alignmentV.6a.sam.C-2S.2  | sort -k1,1 -u > alignmentV.6a.sam.D-all


awk '($22==1)' alignmentV.6a.sam.D-all | awk '{ split($18,a , "S"); print $0 "\t" a[1];}' > alignmentV.6a.sam.Bnext
awk '($26 ~ /M/ )' alignmentV.6a.sam.Bnext | cut -f1-25 > alignmentV.6a.sam.B-MbeforeSnext
awk '($26 !~ /M/ )' alignmentV.6a.sam.Bnext | cut -f1-25 > alignmentV.6a.sam.B-SbeforeMnext
awk '($22==2)'  alignmentV.6a.sam.D-all  > alignmentV.6a.sam.C-2Snext



awk '{print $18}' alignmentV.6a.sam.B-MbeforeSnext | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*S//g' -e 's/[A-Z]//g'  | awk '{tot=0; for(i=1;i<=NF;i++) tot+=$i; print tot}' | awk '{print "0" "\t" $1}' | paste alignmentV.6a.sam.B-MbeforeSnext - > alignmentV.6a.sam.B-MbeforeS.1next
awk '{print $18}' alignmentV.6a.sam.B-MbeforeS.1next | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*I//g' -e 's/[0-9]*M//g' -e 's/[A-Z]//g' | paste alignmentV.6a.sam.B-MbeforeS.1next - > alignmentV.6a.sam.B-MbeforeS.2next

 awk '{print $18}' alignmentV.6a.sam.B-SbeforeMnext | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*I//g' -e 's/[0-9]*M//g' -e 's/[A-Z]//g' | paste alignmentV.6a.sam.B-SbeforeMnext - > alignmentV.6a.sam.B-SbeforeM.1next
 awk '{print $18}' alignmentV.6a.sam.B-SbeforeM.1next | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*S//g' -e 's/[A-Z]//g'  | awk '{tot=0; for(i=1;i<=NF;i++) tot+=$i; print tot}' | awk '{print $1 "\t" "0"}' | paste alignmentV.6a.sam.B-SbeforeM.1next - > alignmentV.6a.sam.B-SbeforeM.2next

 awk '{print $18}'  alignmentV.6a.sam.C-2Snext | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*I//g' -e 's/[0-9]*M//g' -e 's/[A-Z]//g' |  awk '{print $1 }' | paste alignmentV.6a.sam.C-2Snext - > alignmentV.6a.sam.C-2S.1next
awk '{print $18}'  alignmentV.6a.sam.C-2S.1next | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*S//g' -e 's/[A-Z]//g'  | awk '{tot=0; for(i=1;i<=NF;i++) tot+=$i; print tot}' | awk '{print $1 }' | paste alignmentV.6a.sam.C-2S.1next - > alignmentV.6a.sam.C-2S.2anext
 awk '{print $18}'   alignmentV.6a.sam.C-2S.2anext | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*I//g' -e 's/[0-9]*M//g' -e 's/[A-Z]//g' |  awk '{print $2 }' | paste alignmentV.6a.sam.C-2S.2anext - > alignmentV.6a.sam.C-2S.2next

cat alignmentV.6a.sam.B-MbeforeS.2next  alignmentV.6a.sam.B-SbeforeM.2next alignmentV.6a.sam.C-2S.2next | sort -k1,1 -u > alignmentV.6a.sam.D-allnext

cut -f1,26,27,28 alignmentV.6a.sam.D-allnext | paste alignmentV.6a.sam.D-all - | awk '($1==$26)' | cut -f1-25,27- |  sed 's/AS:i://g ' | sed 's/XS:i://g'| awk '{if (($2 ~ /145/) && ($13-$14>=2 && $13-$14<=4)) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" "AS:i:"$13 "\t" "XS:i:""0" "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $19 "\t" $20 "\t" $21 "\t" $22 "\t" $23 "\t" $24 "\t" $25 "\t" $26 "\t" $27 "\t" $28; else print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" "AS:i:"$13 "\t" "XS:i:"$14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $19 "\t" $20 "\t" $21 "\t" $22 "\t" $23 "\t" $24 "\t" $25 "\t" $26 "\t" $27 "\t" $28}'  > alignmentV.7.sam.A

########################
awk '{if ($2==147 || $2==83 || $2==145 || $2==81 || $2==113 || $2==177 ) print $0 "\t",  $2="-"; else print $0 "\t" "+" }'  alignmentV.7.sam.A | awk '{print $0 "\t" $3","$29$4","$6}'  | awk '{print $0"\t"$1":XP:Z:"$30}' | sort -k31,31 > alignmentV.7.sam.B

if [[ $mode == *sensitive* ]]; then
	awk '($2>=255 || $2=="131" || $2=="0" || $2=="16")' alignmentV.3.sam.ext | awk '($15 ~ /XP:Z:/)' | awk '{split ($15, a, ";"); print ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, a[1], $1":"a[1])}' OFS='\t' | sort -k16,16 -u  > alignmentV.3.samA
else
	awk '($2>=255 || $2=="131" || $2=="0" || $2=="16")' alignmentV.3.sam | awk '($15 ~ /XP:Z:/)' | awk '{split ($15, a, ";"); print ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, a[1], $1":"a[1])}' OFS='\t' | sort -k16,16 -u  > alignmentV.3.samA
fi
awk '{print $31 }' alignmentV.7.sam.B | grep -Fwf - alignmentV.3.samA | sort -k1,1 >  alignmentV.3.sam.temp1

sort -k1,1 alignmentV.7.sam.B | cut -f1-31 | paste - alignmentV.3.sam.temp1 | cut -f1-46 |awk '{if ($1==$32) print $0; else print "ERROR!!! ANALYSIS FAILED."}' > alignmentV.8.sam.A
############
awk '{print $1 "\t" $3 "\t" $4 "\t" $29 "\t" $6 "\t" $10 "\t" $12 "\t" $14 "\t" $34 "\t" $35 "\t" $15 "\t" $37 "\t" $32 "\t" $43 "\t" $45 "\t" $13 "\t" $44 "\t" $23 "\t" $24 "\t" $25 "\t" $26 "\t" $27 "\t" $28}' alignmentV.8.sam.A  | awk '{if($11 ~ /+/) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" "+" "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $19 "\t" $20 "\t" $21 "\t" $22 "\t" $23; else print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" "-" "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $19 "\t" $20 "\t" $21 "\t" $22 "\t" $23 }' > alignmentV.8.sam.B

awk '{if ($18!=0 && $20!=0 && ($21==0 || $23==0)) print $0 }' alignmentV.8.sam.B |awk '{ split($12,a , "S"); print $0 "\t" a[1];}'| awk '{if (($24 !~ /M/) && (($4=="+" && $11=="+") || ($4=="-" && $11=="-"))) {OFS= "\t"; $18=0; print $0}; if (($24 !~ /M/) && (($4=="+" && $11=="-") || ($4=="-" && $11=="+"))) {OFS= "\t"; $20=0; print $0}; if (($24 ~ /M/) && (($4=="+" && $11=="+") || ($4=="-" && $11=="-"))) {OFS= "\t"; $20=0; print $0}; if (($24 ~ /M/) && (($4=="+" && $11=="-") || ($4=="-" && $11=="+"))) {OFS= "\t"; $18=0; print $0};}' | cut -f1-23 > alignmentV.9.sam.A

awk '{if ($21!=0 && $23!=0 && ($18==0 || $20==0)) print $0 }' alignmentV.8.sam.B | awk '{print $1 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $6 "\t" $14 "\t" $15 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $13 "\t" $7 "\t" $8 "\t" $17 "\t" $16 "\t" $21 "\t" $22 "\t" $23 "\t" $18 "\t" $19 "\t" $20}' > alignmentV.9.sam.B0

awk '{ split($12,a , "S"); print $0 "\t" a[1];}' alignmentV.9.sam.B0 | awk '{if (($24 !~ /M/) && (($4=="+" && $11=="+") || ($4=="-" && $11=="-"))) {OFS= "\t"; $18=0; print $0}; if (($24 !~ /M/) && (($4=="+" && $11=="-") || ($4=="-" && $11=="+"))) {OFS= "\t"; $20=0; print $0}; if (($24 ~ /M/) && (($4=="+" && $11=="+") || ($4=="-" && $11=="-"))) {OFS= "\t"; $20=0; print $0}; if (($24 ~ /M/) && (($4=="+" && $11=="-") || ($4=="-" && $11=="+"))) {OFS= "\t"; $18=0; print $0};}' | cut -f1-23 > alignmentV.9.sam.B


awk '{if (($18==0 || $20==0) && ($21==0 || $23==0)) print $0 }' alignmentV.8.sam.B > alignmentV.9.sam.D

if [[ $mode == *sensitive* ]]; then
	#if double S in both reads then link it to SENS.
	awk '{if ($18!=0 && $20!=0 && $21!=0 && $23!=0) print $0 }' alignmentV.8.sam.B  | awk '(($18<=10 || $20<=10) && ($21<=10 || $23<=10))' | awk '{if ($18>=11 && $20<=10) {OFS= "\t"; $20=0; print $0}; if ($20>=11 && $18<=10) {OFS= "\t"; $18=0; print $0}}' |  awk '{if ($23>=11 && $21<=10) {OFS= "\t"; $21=0; print $0}; if ($21>=11 && $23<=10) {OFS= "\t"; $23=0; print $0}}' > alignmentV.9.sam.C
	cat alignmentV.9.sam.A alignmentV.9.sam.B alignmentV.9.sam.C alignmentV.9.sam.D > alignmentV.10.sam.A
fi
cat alignmentV.9.sam.A alignmentV.9.sam.B alignmentV.9.sam.D > alignmentV.10.sam.A

###########################
echo " "
echo "Fusion positions processing and parameter based estimation and filtering..."

awk '{if ($18==0 && $21==0) print $0 "\t" $22-$20; if ($18==0 && $23==0) print $0 "\t" $22-$20; if ($20==0 && $21==0) print $0 "\t" $22-$18; if ($20==0 && $23==0) print $0 "\t" $22-$18;}' alignmentV.10.sam.A  > alignmentV.11.sam.A


 awk '{if (($18==0 && $21==0) && ($4=="+")) print $0 "\t" $3+$19-1; if (($18==0 && $21==0) && ($4=="-"))  print $0 "\t" $3+$19-1 }' alignmentV.11.sam.A | awk '{if ($11=="+") print $0 "\t" $10+$22-1; if ($11=="-")  print $0 "\t" $10+$22-1}' > alignmentV.12.sam.A

awk '{ print $1 "\t" $2 "\t" $25 "\t" $4 "\t" $9 "\t" $11 "\t" $26 "\t" $6 "\t" $19 "\t" $22 "\t" $24 "\t" $7 "\t" $14 "\t" $16 "\t" $17 "\t" $8 "\t" $15 "\t" $5}' alignmentV.12.sam.A | awk '{if ($11 <= 0) print $0 "\t" "0" "\t" $11 ; if ($11>=1) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7-$11 "\t" $8 "\t" $9 "\t" $10-$11 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $11 "\t" "0"}'  | awk -F"\t" '{gsub(/-/,"",$20)}1' OFS="\t" >  alignmentV.14.sam.1

awk '{if (($18==0 && $23==0) && ($4=="+")) print $0 "\t" $3+$19-1; if (($18==0 && $23==0) && ($4=="-"))  print $0 "\t" $3+$19-1 }' alignmentV.11.sam.A | awk '{if ($11=="+") print $0 "\t" $10; if ($11=="-")  print $0 "\t" $10 }' > alignmentV.12.sam.B

awk '{ print $1 "\t" $2 "\t" $25 "\t" $4 "\t" $9 "\t" $11 "\t" $26 "\t" $6 "\t" $19 "\t" $22 "\t" $24 "\t" $7 "\t" $14 "\t" $16 "\t" $17 "\t" $8 "\t" $15 "\t" $5}' alignmentV.12.sam.B | awk '{if ($11 <= 0) print $0 "\t" "0" "\t" $11 ; if ($11>=1) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7+$11 "\t" $8 "\t" $9 "\t" $10-$11 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $11 "\t" "0"}'  | awk -F"\t" '{gsub(/-/,"",$20)}1' OFS="\t" >  alignmentV.14.sam.2

awk '{if (($20==0 && $23==0) && ($4=="+")) print $0 "\t" $3; if (($20==0 && $23==0) && ($4=="-"))  print $0 "\t" $3 }' alignmentV.11.sam.A | awk '{if ($11=="+") print $0 "\t" $10; if ($11=="-")  print $0 "\t" $10 }' > alignmentV.12.sam.C

awk '{ print $1 "\t" $2 "\t" $25 "\t" $4 "\t" $9 "\t" $11 "\t" $26 "\t" $6 "\t" $19 "\t" $22 "\t" $24 "\t" $7 "\t" $14 "\t" $16 "\t" $17 "\t" $8 "\t" $15 "\t" $5}' alignmentV.12.sam.C | awk '{if ($11 <= 0) print $0 "\t" "0" "\t" $11 ; if ($11>=1) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7+$11 "\t" $8 "\t" $9 "\t" $10-$11 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $11 "\t" "0"}'  | awk -F"\t" '{gsub(/-/,"",$20)}1' OFS="\t" >  alignmentV.14.sam.3

awk '{if (($20==0 && $21==0) && ($4=="+")) print $0 "\t" $3; if (($20==0 && $21==0) && ($4=="-"))  print $0 "\t" $3 }' alignmentV.11.sam.A | awk '{if ($11=="+") print $0 "\t" $10+$22-1; if ($11=="-")  print $0 "\t" $10+$22-1 }' > alignmentV.12.sam.D

awk '{ print $1 "\t" $2 "\t" $25 "\t" $4 "\t" $9 "\t" $11 "\t" $26 "\t" $6 "\t" $19 "\t" $22 "\t" $24 "\t" $7 "\t" $14 "\t" $16 "\t" $17 "\t" $8 "\t" $15 "\t" $5}' alignmentV.12.sam.D | awk '{if ($11 <= 0) print $0 "\t" "0" "\t" $11 ; if ($11>=1) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7-$11 "\t" $8 "\t" $9 "\t" $10-$11 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $11 "\t" "0"}'  | awk -F"\t" '{gsub(/-/,"",$20)}1' OFS="\t" >  alignmentV.14.sam.4
cat alignmentV.14.sam.1 alignmentV.14.sam.2 alignmentV.14.sam.3 alignmentV.14.sam.4 > alignmentV.14.sam.A


awk -v distVecVec=$distVecVec '($20<=distVecVec)' alignmentV.14.sam.A  > alignmentV.15.sam.A


awk -v opVecVec=$opVecVec '($19<=opVecVec)' alignmentV.15.sam.A > alignmentV.16.sam.A


awk -v minMapSpanVec=$minMapSpanVec '($9>=minMapSpanVec && $10>=minMapSpanVec)' alignmentV.16.sam.A > alignmentV.17.sam.A

awk -F"\t" '{gsub(/AS:i:/,"",$0)}1' OFS="\t" alignmentV.17.sam.A | awk -F"\t" '{gsub(/XS:i:/,"",$0)}1' OFS="\t" | awk '{print $0 "\t" $16/$14 "\t" $17/$15}' | cut -f1-10,12,13,18,19,20,21,22 > alignmentV.18.sam.A
awk -v UMthresholdVec=$UMthresholdVec '{if ($16>=UMthresholdVec) print $0 "\t" "Non-unique"; else print $0 "\t" "Unique"}' alignmentV.18.sam.A | awk -v UMthresholdVec=$UMthresholdVec '{if ($17>=UMthresholdVec) print $0 "\t" "Non-unique"; else print $0 "\t" "Unique"}' > alignmentV.18.sam.B

awk -F"\t" '{gsub(/NM:i:/,"",$0)}1' OFS="\t" alignmentV.18.sam.B |awk '{print $0 "\t" ((($9-$11)/$9)*100) "\t" (((($10+$14)-$12)/($10+$14))*100)}' | cut -f1-10,14- | awk -v idenVecVec=$idenVecVec '($17>=idenVecVec && $18>=idenVecvec)' > alignmentV.19.sam.A

awk '{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $7 "\t" $6 "\t" $1 "\t "$8 "\t" $9 "\t" $10 "\t" $13 "\t" $15 "\t" $14 "\t" $16 "\t" $17 "\t" $18 "\t" $11 "\t" $12}' alignmentV.19.sam.A | awk '($11<=1 && $13<=1)' > alignmentV.20.sam.A 
######################
echo " "
echo "Duplicate reads removal in process..."

cut -f7 $outDir/alignmentV.20.sam.A > $outDir/alignmentV.21.sam.temp.ids
awk 'NR==FNR{tgts[$1]; next} $1 in tgts'  $outDir/alignmentV.21.sam.temp.ids $outDir/alignmentV.3.sam > $outDir/alignmentV.21.sam.temp.ids.full
$samtools view -bT $vecRef $outDir/alignmentV.21.sam.temp.ids.full  > $outDir/alignmentV.21.sam.temp.ids.full.bam
$samtools sort $outDir/alignmentV.21.sam.temp.ids.full.bam -o  $outDir/alignmentV.21.sam.temp.ids.full.sorted.bam
$samtools rmdup $outDir/alignmentV.21.sam.temp.ids.full.sorted.bam $outDir/alignmentV.21.sam.temp.ids.full.sorted.nodup.bam
$samtools view  $outDir/alignmentV.21.sam.temp.ids.full.sorted.nodup.bam  |cut -f1 |sort | uniq > $outDir/alignmentV.21.sam.temp.ids.full.sorted.nodup.ids


awk 'NR==FNR{tgts[$1]; next} $7 in tgts' $outDir/alignmentV.21.sam.temp.ids.full.sorted.nodup.ids $outDir/alignmentV.20.sam.A > $outDir/alignmentV.22.sam.A

echo " "
echo "Final breakpoints processing per viral vector..."
echo " "
echo "Generating final result files for vector-vector fusions..."
mv alignmentV.22.sam.A alignmentV.23.sam.A

#split file into multiple based on vector
if [ -s "$outDir/alignmentV.23.sam.A" ]; then
rm -f *.IntraVector_BreakPoints*
sort -k4,4 -k1,1 -k2,2n $outDir/alignmentV.23.sam.A | cut -f1-10,12,14- | awk '{print>$4".IntraVector_BreakPoints"}'
sed -i 's/\t/,/g' *.IntraVector_BreakPoints
echo "VectorName,BreakPointPosition1,Strand1,VectorName,BreakPointPosition2,Strand2,ReadID,Sequence,Span1,Span2,Feature1,Feature2,IdentityPerc1,IdentityPerc2,OverlapFusion,DistanceFusion" > $outDir/header.fusion
for f in *IntraVector_BreakPoints ; do cat header.fusion $f >>$f.csv; done

#for f in *IntraVector_BreakPoints; do awk -F ','  '{print $1 "\t" $2 "\n" $4 "\t" $5}' $f |  sort | uniq -c |  sed -e 's/^[ \t]*//' | awk -F ' ' '{print $2 "\t" $3 "\t" $3-1 "\t" "FF" "\t" $1}' | sort -k2,2n -k3,3n >> $f.IGV ; done

mv alignmentV.21.sam.temp.ids  vecToRemoveIDS
#mv alignmentV.3.sam  alignmentV_vecVecFusion.sam
else 
echo " "
echo "NO VECTOR-VECTOR FUSION READS FOUND."
fi
echo " "
echo "    Vector-Vector Fusion Analysis Finished"
echo "..........................................."
rm -f file.log   *[*]* alignmentV.* header*  *.IntraVector_BreakPoints *Hit.split.txt.cluster* alignmentV_vecVecFusion.sam
echo "Vector-Vector Fusion Analysis finished" >> log
date >> log
############################################
