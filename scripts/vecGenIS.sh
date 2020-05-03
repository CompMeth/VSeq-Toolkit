outDir=$1
file1=$2
file2=$3
vecRef=$4
stringencyVecGen=$5
UMthresholdVecGen=$6
minMapSpanVecGen=$7
aligner=$8
samtools=$9
vecGenRef=${10}
distVecGen=${11}
opVecGen=${12}
idenVecGen=${13}
clusterRange=${14}
bedtools=${15}
annoTable=${16}
bin=${17}
vecVecFusionCheck=${18}
mode=${19}
contCheck=${20}

cd $outDir
echo "Host-Vector Fusion Analysis started" >> log
date >> log

echo " "
if [[ $contCheck == *true* && $vecVecFusionCheck == *true* ]]; then
        echo "Alignment with viral vector and genome reference in process..."
        $aligner mem -M -t 10 -T 15 $vecGenRef $outDir/ContAnalysis_WithoutContaminantsReadPairs.R1.fastq.gz  $outDir/ContAnalysis_WithoutContaminantsReadPairs.R2.fastq.gz > $outDir/alignment.2.sam.basic0
        grep -v -Fwf vecToRemoveIDS $outDir/alignment.2.sam.basic0 > $outDir/alignment.2.sam.basic
fi

if [[ $contCheck == *false* && $vecVecFusionCheck == *true* ]]; then
	echo "Alignment with viral vector and genome reference in process..."
	$aligner mem -M -t 10 -T 15 $vecGenRef $outDir/*pair1.fastq  $outDir/*pair2.fastq > $outDir/alignment.2.sam.basic0	
	grep -v -Fwf vecToRemoveIDS $outDir/alignment.2.sam.basic0 > $outDir/alignment.2.sam.basic
fi

if [[ $contCheck == *false* && $vecVecFusionCheck == *false* ]];then
		echo ""
		echo "Alignment with viral vector and genome reference in process..."
		$aligner mem -M -t 10 -T 15 $vecGenRef $outDir/*pair1.fastq  $outDir/*pair2.fastq > $outDir/alignment.2.sam.basic
fi

awk '($6 ~ /S/)' $outDir/alignment.2.sam.basic > $outDir/alignment.2.sam


if [[ $mode == *sensitive* ]]; then
	echo ""
        echo "Sensitive mode selected"
	echo "Processing started..."
        mv alignment.2.sam alignmentSens.sam
        bash $bin/sensModule.sh $outDir $vecRef $aligner
        mv alignmentSens.sam.1 alignment.2.sam
        mv alignmentSens.3.sam.ext alignment.3.sam.ext
fi

mv $outDir/alignment.2.sam $outDir/alignment.3.sam

echo " "
echo "Selection of potential fusion reads in process..."
echo " "
if [[ $stringencyVecGen == *high* ]]; then
        awk '($2=="99" || $2=="147" || $2=="83" || $2=="163")' $outDir/alignment.3.sam  | grep 'XP:Z:' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 | sed 's/XP:Z://g' > $outDir/alignment.4.sam
fi
if [[ $stringencyVecGen == *medium* ]]; then
        awk '($2=="99" || $2=="147" || $2=="83" || $2=="163" || $2=="81" || $2=="161" || $2=="97" || $2=="145")' $outDir/alignment.3.sam  | grep 'XP:Z:' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 | sed 's/XP:Z://g' > $outDir/alignment.4.sam
fi
if [[ $stringencyVecGen == *low* ]]; then
	awk '($2=="99" || $2=="147" || $2=="83" || $2=="163" || $2=="81" || $2=="161" || $2=="97" || $2=="145"  || $2=="65"  || $2=="129" || $2=="113" || $2=="177")' $outDir/alignment.3.sam  | grep 'XP:Z:' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 | sed 's/XP:Z://g' > $outDir/alignment.4.sam
fi
if [[ $stringencyVecGen == *null* ]]; then
        awk '($2=="99" || $2=="147" || $2=="83" || $2=="163" || $2=="81" || $2=="161" || $2=="97" || $2=="145"  || $2=="65"  || $2=="129" || $2=="113" || $2=="177" || $2=="131" || $2=="387" )' $outDir/alignment.3.sam  | grep 'XP:Z:' | cut -f1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 | sed 's/XP:Z://g' > $outDir/alignment.4.sam
fi

echo " "
echo "Post-processing of extracted reads..."
echo " "
cat $outDir/alignment.4.sam  |  cut -f15 | cut -d ';' -f1 | paste $outDir/alignment.4.sam - | cut -f1-14,16  | awk '{print $0";"}' > $outDir/alignment.4.sam.temp00
cut -f15  $outDir/alignment.4.sam.temp00|  awk '{print gsub (/;/, "")}' | paste $outDir/alignment.4.sam.temp00 - | awk '($16>=2)' | cut -f1 | grep -v -Fwf - $outDir/alignment.4.sam.temp00 > $outDir/alignment.4.sam.temp1
mv $outDir/alignment.4.sam.temp00 $outDir/alignment.4.sam

cut -f15 $outDir/alignment.4.sam  | sed 's/,/\t/g' > $outDir/alignment.5.sam
paste $outDir/alignment.4.sam $outDir/alignment.5.sam | awk '($3!=$16)' > $outDir/alignment.6.sam
awk '($3 ~ /chr/ || $16 ~ /chr/)' $outDir/alignment.6.sam  | awk '($3 !~ /chr/ || $16 !~ /chr/)' > $outDir/alignment.6a.sam

################

awk '{print $0, gsub(/S/,"\t",$6) }' alignment.6a.sam | awk '{print $0, gsub(/S/,"\t",$18) }' | sed 's/ /\t/g' > alignment.6a.sam.A


awk '($21==1)' alignment.6a.sam.A | awk '{ split($6,a , "S"); print $0 "\t" a[1];}' > alignment.6a.sam.B
awk '($23 ~ /M/ )' alignment.6a.sam.B | cut -f1-22 > alignment.6a.sam.B-MbeforeS
awk '($23 !~ /M/ )' alignment.6a.sam.B | cut -f1-22 > alignment.6a.sam.B-SbeforeM

awk '($21==2)' alignment.6a.sam.A  > alignment.6a.sam.C-2S
########################


awk '{print $6}' alignment.6a.sam.B-MbeforeS | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*S//g' -e 's/[A-Z]//g'  | awk '{tot=0; for(i=1;i<=NF;i++) tot+=$i; print tot}' | awk '{print "0" "\t" $1}' | paste alignment.6a.sam.B-MbeforeS - > alignment.6a.sam.B-MbeforeS.1
awk '{print $6}' alignment.6a.sam.B-MbeforeS.1 | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*I//g' -e 's/[0-9]*M//g' -e 's/[A-Z]//g' | paste alignment.6a.sam.B-MbeforeS.1 - > alignment.6a.sam.B-MbeforeS.2


 awk '{print $6}' alignment.6a.sam.B-SbeforeM | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*I//g' -e 's/[0-9]*M//g' -e 's/[A-Z]//g' | paste alignment.6a.sam.B-SbeforeM - > alignment.6a.sam.B-SbeforeM.1
 awk '{print $6}' alignment.6a.sam.B-SbeforeM.1 | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*S//g' -e 's/[A-Z]//g'  | awk '{tot=0; for(i=1;i<=NF;i++) tot+=$i; print tot}' | awk '{print $1 "\t" "0"}' | paste alignment.6a.sam.B-SbeforeM.1 - > alignment.6a.sam.B-SbeforeM.2 


 awk '{print $6}' alignment.6a.sam.C-2S | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*I//g' -e 's/[0-9]*M//g' -e 's/[A-Z]//g' |  awk '{print $1 }' | paste alignment.6a.sam.C-2S - > alignment.6a.sam.C-2S.1
awk '{print $6}' alignment.6a.sam.C-2S.1 | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*S//g' -e 's/[A-Z]//g'  | awk '{tot=0; for(i=1;i<=NF;i++) tot+=$i; print tot}' | awk '{print $1 }' | paste alignment.6a.sam.C-2S.1 - > alignment.6a.sam.C-2S.2a
 awk '{print $6}'  alignment.6a.sam.C-2S.2a | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*I//g' -e 's/[0-9]*M//g' -e 's/[A-Z]//g' |  awk '{print $2 }' | paste alignment.6a.sam.C-2S.2a - > alignment.6a.sam.C-2S.2

cat alignment.6a.sam.B-MbeforeS.2  alignment.6a.sam.B-SbeforeM.2 alignment.6a.sam.C-2S.2  | sort -k1,1 -u > alignment.6a.sam.D-all
########################

awk '($22==1)' alignment.6a.sam.D-all | awk '{ split($18,a , "S"); print $0 "\t" a[1];}' > alignment.6a.sam.Bnext
awk '($26 ~ /M/ )' alignment.6a.sam.Bnext | cut -f1-25 > alignment.6a.sam.B-MbeforeSnext
awk '($26 !~ /M/ )' alignment.6a.sam.Bnext | cut -f1-25 > alignment.6a.sam.B-SbeforeMnext

awk '($22==2)'  alignment.6a.sam.D-all  > alignment.6a.sam.C-2Snext
########################


awk '{print $18}' alignment.6a.sam.B-MbeforeSnext | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*S//g' -e 's/[A-Z]//g'  | awk '{tot=0; for(i=1;i<=NF;i++) tot+=$i; print tot}' | awk '{print "0" "\t" $1}' | paste alignment.6a.sam.B-MbeforeSnext - > alignment.6a.sam.B-MbeforeS.1next
awk '{print $18}' alignment.6a.sam.B-MbeforeS.1next | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*I//g' -e 's/[0-9]*M//g' -e 's/[A-Z]//g' | paste alignment.6a.sam.B-MbeforeS.1next - > alignment.6a.sam.B-MbeforeS.2next


 awk '{print $18}' alignment.6a.sam.B-SbeforeMnext | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*I//g' -e 's/[0-9]*M//g' -e 's/[A-Z]//g' | paste alignment.6a.sam.B-SbeforeMnext - > alignment.6a.sam.B-SbeforeM.1next
 awk '{print $18}' alignment.6a.sam.B-SbeforeM.1next | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*S//g' -e 's/[A-Z]//g'  | awk '{tot=0; for(i=1;i<=NF;i++) tot+=$i; print tot}' | awk '{print $1 "\t" "0"}' | paste alignment.6a.sam.B-SbeforeM.1next - > alignment.6a.sam.B-SbeforeM.2next


 awk '{print $18}'  alignment.6a.sam.C-2Snext | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*I//g' -e 's/[0-9]*M//g' -e 's/[A-Z]//g' |  awk '{print $1 }' | paste alignment.6a.sam.C-2Snext - > alignment.6a.sam.C-2S.1next
awk '{print $18}'  alignment.6a.sam.C-2S.1next | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*S//g' -e 's/[A-Z]//g'  | awk '{tot=0; for(i=1;i<=NF;i++) tot+=$i; print tot}' | awk '{print $1 }' | paste alignment.6a.sam.C-2S.1next - > alignment.6a.sam.C-2S.2anext
 awk '{print $18}'   alignment.6a.sam.C-2S.2anext | sed -e 's/[A-Z]/& /g' -e  's/[0-9]*D//g' -e 's/[0-9]*I//g' -e 's/[0-9]*M//g' -e 's/[A-Z]//g' |  awk '{print $2 }' | paste alignment.6a.sam.C-2S.2anext - > alignment.6a.sam.C-2S.2next

cat alignment.6a.sam.B-MbeforeS.2next  alignment.6a.sam.B-SbeforeM.2next alignment.6a.sam.C-2S.2next | sort -k1,1 -u > alignment.6a.sam.D-allnext

cut -f1,26,27,28 alignment.6a.sam.D-allnext | paste alignment.6a.sam.D-all - | awk '($1==$26)' | cut -f1-25,27-  |  sed 's/AS:i://g ' | sed 's/XS:i://g'| awk '{if (($2 ~ /145/) && ($13-$14>=2 && $13-$14<=4)) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" "AS:i:"$13 "\t" "XS:i:""0" "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $19 "\t" $20 "\t" $21 "\t" $22 "\t" $23 "\t" $24 "\t" $25 "\t" $26 "\t" $27 "\t" $28; else print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" "AS:i:"$13 "\t" "XS:i:"$14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $19 "\t" $20 "\t" $21 "\t" $22 "\t" $23 "\t" $24 "\t" $25 "\t" $26 "\t" $27 "\t" $28}' > alignment.7.sam.A
########################
awk '{if ($2==147 || $2==83 || $2==145 || $2==81 || $2==113 || $2==177 ) print $0 "\t",  $2="-"; else print $0 "\t" "+" }'  alignment.7.sam.A | awk '{print $0 "\t" $3","$29$4","$6}'  | awk '{print $0"\t"$1":XP:Z:"$30}' | sort -k31,31 > alignment.7.sam.B

if [[ $mode == *sensitive* ]]; then
	awk '($2>=255 || $2=="131" || $2=="0" || $2=="16")' alignment.3.sam.ext | awk '($15 ~ /XP:Z:/)' | awk '{split ($15, a, ";"); print ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, a[1], $1":"a[1])}' OFS='\t' | sort -k16,16 -u  > alignment.3.samA
else
	awk '($2>=255 || $2=="131" || $2=="0" || $2=="16")' alignment.3.sam | awk '($15 ~ /XP:Z:/)' | awk '{split ($15, a, ";"); print ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, a[1], $1":"a[1])}' OFS='\t' | sort -k16,16 -u  > alignment.3.samA
fi
awk '{print $31 }' alignment.7.sam.B | grep -Fwf - alignment.3.samA | sort -k1,1 >  alignment.3.sam.temp1

sort -k1,1 alignment.7.sam.B | cut -f1-31 | paste - alignment.3.sam.temp1 | cut -f1-46 |awk '{if ($1==$32) print $0; else print "ERROR!!! ANALYSIS FAILED."}' > alignment.8.sam.A
########################

 awk '{print $1 "\t" $3 "\t" $4 "\t" $29 "\t" $6 "\t" $10 "\t" $12 "\t" $14 "\t" $34 "\t" $35 "\t" $15 "\t" $37 "\t" $32 "\t" $43 "\t" $45 "\t" $13 "\t" $44 "\t" $23 "\t" $24 "\t" $25 "\t" $26 "\t" $27 "\t" $28}' alignment.8.sam.A  | awk '{if($11 ~ /+/) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" "+" "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $19 "\t" $20 "\t" $21 "\t" $22 "\t" $23; else print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" "-" "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $19 "\t" $20 "\t" $21 "\t" $22 "\t" $23 }' > alignment.8.sam.B

awk '{if ($18!=0 && $20!=0 && ($21==0 || $23==0)) print $0 }' alignment.8.sam.B |awk '{ split($12,a , "S"); print $0 "\t" a[1];}'| awk '{if (($24 !~ /M/) && (($4=="+" && $11=="+") || ($4=="-" && $11=="-"))) {OFS= "\t"; $18=0; print $0}; if (($24 !~ /M/) && (($4=="+" && $11=="-") || ($4=="-" && $11=="+"))) {OFS= "\t"; $20=0; print $0}; if (($24 ~ /M/) && (($4=="+" && $11=="+") || ($4=="-" && $11=="-"))) {OFS= "\t"; $20=0; print $0}; if (($24 ~ /M/) && (($4=="+" && $11=="-") || ($4=="-" && $11=="+"))) {OFS= "\t"; $18=0; print $0};}' | cut -f1-23 > alignment.9.sam.A
#if 2nd cigar has 2 S then this;
awk '{if ($21!=0 && $23!=0 && ($18==0 || $20==0)) print $0 }' alignment.8.sam.B | awk '{print $1 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $6 "\t" $14 "\t" $15 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $13 "\t" $7 "\t" $8 "\t" $17 "\t" $16 "\t" $21 "\t" $22 "\t" $23 "\t" $18 "\t" $19 "\t" $20}' > alignment.9.sam.B0

awk '{ split($12,a , "S"); print $0 "\t" a[1];}' alignment.9.sam.B0 | awk '{if (($24 !~ /M/) && (($4=="+" && $11=="+") || ($4=="-" && $11=="-"))) {OFS= "\t"; $18=0; print $0}; if (($24 !~ /M/) && (($4=="+" && $11=="-") || ($4=="-" && $11=="+"))) {OFS= "\t"; $20=0; print $0}; if (($24 ~ /M/) && (($4=="+" && $11=="+") || ($4=="-" && $11=="-"))) {OFS= "\t"; $20=0; print $0}; if (($24 ~ /M/) && (($4=="+" && $11=="-") || ($4=="-" && $11=="+"))) {OFS= "\t"; $18=0; print $0};}' | cut -f1-23 > alignment.9.sam.B


awk '{if (($18==0 || $20==0) && ($21==0 || $23==0)) print $0 }' alignment.8.sam.B > alignment.9.sam.D

if [[ $mode == *sensitive* ]]; then
        #if double S in both reads then link it to SENS.
        awk '{if ($18!=0 && $20!=0 && $21!=0 && $23!=0) print $0 }' alignment.8.sam.B  | awk '(($18<=10 || $20<=10) && ($21<=10 || $23<=10))' | awk '{if ($18>=11 && $20<=10) {OFS= "\t"; $20=0; print $0}; if ($20>=11 && $18<=10) {OFS= "\t"; $18=0; print $0}}' |  awk '{if ($23>=11 && $21<=10) {OFS= "\t"; $21=0; print $0}; if ($21>=11 && $23<=10) {OFS= "\t"; $23=0; print $0}}' > alignment.9.sam.C       
        cat alignment.9.sam.A alignment.9.sam.B alignment.9.sam.C alignment.9.sam.D > alignment.10.sam.A
fi
cat alignment.9.sam.A alignment.9.sam.B alignment.9.sam.D > alignment.10.sam.A
###########################
echo " "
echo "Fusion positions processing and parameter based estimation and filtering..."
#calculate overlap or distilled
awk '{if ($18==0 && $21==0) print $0 "\t" $22-$20; if ($18==0 && $23==0) print $0 "\t" $22-$20; if ($20==0 && $21==0) print $0 "\t" $22-$18; if ($20==0 && $23==0) print $0 "\t" $22-$18;}' alignment.10.sam.A  > alignment.11.sam.A


 awk '{if (($18==0 && $21==0) && ($4=="+")) print $0 "\t" $3+$19-1; if (($18==0 && $21==0) && ($4=="-"))  print $0 "\t" $3+$19-1 }' alignment.11.sam.A | awk '{if ($11=="+") print $0 "\t" $10+$22-1; if ($11=="-")  print $0 "\t" $10+$22-1}' > alignment.12.sam.A

awk '{ print $1 "\t" $2 "\t" $25 "\t" $4 "\t" $9 "\t" $11 "\t" $26 "\t" $6 "\t" $19 "\t" $22 "\t" $24 "\t" $7 "\t" $14 "\t" $16 "\t" $17 "\t" $8 "\t" $15 "\t" $5}' alignment.12.sam.A | awk '{if ($11 <= 0) print $0 "\t" "0" "\t" $11 ; if ($11>=1) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7-$11 "\t" $8 "\t" $9 "\t" $10-$11 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $11 "\t" "0"}'  | awk -F"\t" '{gsub(/-/,"",$20)}1' OFS="\t" >  alignment.14.sam.1

awk '{if (($18==0 && $23==0) && ($4=="+")) print $0 "\t" $3+$19-1; if (($18==0 && $23==0) && ($4=="-"))  print $0 "\t" $3+$19-1 }' alignment.11.sam.A | awk '{if ($11=="+") print $0 "\t" $10; if ($11=="-")  print $0 "\t" $10 }' > alignment.12.sam.B

awk '{ print $1 "\t" $2 "\t" $25 "\t" $4 "\t" $9 "\t" $11 "\t" $26 "\t" $6 "\t" $19 "\t" $22 "\t" $24 "\t" $7 "\t" $14 "\t" $16 "\t" $17 "\t" $8 "\t" $15 "\t" $5}' alignment.12.sam.B | awk '{if ($11 <= 0) print $0 "\t" "0" "\t" $11 ; if ($11>=1) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7+$11 "\t" $8 "\t" $9 "\t" $10-$11 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $11 "\t" "0"}'  | awk -F"\t" '{gsub(/-/,"",$20)}1' OFS="\t" >  alignment.14.sam.2

awk '{if (($20==0 && $23==0) && ($4=="+")) print $0 "\t" $3; if (($20==0 && $23==0) && ($4=="-"))  print $0 "\t" $3 }' alignment.11.sam.A | awk '{if ($11=="+") print $0 "\t" $10; if ($11=="-")  print $0 "\t" $10 }' > alignment.12.sam.C

awk '{ print $1 "\t" $2 "\t" $25 "\t" $4 "\t" $9 "\t" $11 "\t" $26 "\t" $6 "\t" $19 "\t" $22 "\t" $24 "\t" $7 "\t" $14 "\t" $16 "\t" $17 "\t" $8 "\t" $15 "\t" $5}' alignment.12.sam.C | awk '{if ($11 <= 0) print $0 "\t" "0" "\t" $11 ; if ($11>=1) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7+$11 "\t" $8 "\t" $9 "\t" $10-$11 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $11 "\t" "0"}'  | awk -F"\t" '{gsub(/-/,"",$20)}1' OFS="\t" >  alignment.14.sam.3

awk '{if (($20==0 && $21==0) && ($4=="+")) print $0 "\t" $3; if (($20==0 && $21==0) && ($4=="-"))  print $0 "\t" $3 }' alignment.11.sam.A | awk '{if ($11=="+") print $0 "\t" $10+$22-1; if ($11=="-")  print $0 "\t" $10+$22-1 }' > alignment.12.sam.D

awk '{ print $1 "\t" $2 "\t" $25 "\t" $4 "\t" $9 "\t" $11 "\t" $26 "\t" $6 "\t" $19 "\t" $22 "\t" $24 "\t" $7 "\t" $14 "\t" $16 "\t" $17 "\t" $8 "\t" $15 "\t" $5}' alignment.12.sam.D | awk '{if ($11 <= 0) print $0 "\t" "0" "\t" $11 ; if ($11>=1) print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7-$11 "\t" $8 "\t" $9 "\t" $10-$11 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $11 "\t" "0"}'  | awk -F"\t" '{gsub(/-/,"",$20)}1' OFS="\t" >  alignment.14.sam.4

cat alignment.14.sam.1 alignment.14.sam.2 alignment.14.sam.3 alignment.14.sam.4 > alignment.14.sam.A


awk -v distVecGen=$distVecGen '($20<=distVecGen)' alignment.14.sam.A  > alignment.15.sam.A


awk -v opVecGen=$opVecGen '($19<=opVecGen)' alignment.15.sam.A > alignment.16.sam.A


awk -v minMapSpanVecGen=$minMapSpanVecGen '($9>=minMapSpanVecGen && $10>=minMapSpanVecGen)' alignment.16.sam.A > alignment.17.sam.A


awk -F"\t" '{gsub(/AS:i:/,"",$0)}1' OFS="\t" alignment.17.sam.A | awk -F"\t" '{gsub(/XS:i:/,"",$0)}1' OFS="\t" | awk '{print $0 "\t" $16/$14 "\t" $17/$15}' | cut -f1-10,12,13,18,19,20,21,22 > alignment.18.sam.A
 awk -v UMthresholdVecGen=$UMthresholdVecGen '{if ($16>=UMthresholdVecGen) print $0 "\t" "Non-uniqueHit"; else print $0 "\t" "UniqueHit"}' alignment.18.sam.A | awk -v UMthresholdVecGen=$UMthresholdVecGen '{if ($17>=UMthresholdVecGen) print $0 "\t" "Non-uniqueHit"; else print $0 "\t" "UniqueHit"}' > alignment.18.sam.B


 awk -F"\t" '{gsub(/NM:i:/,"",$0)}1' OFS="\t" alignment.18.sam.B |awk '{print $0 "\t" ((($9-$11)/$9)*100) "\t" (((($10+$14)-$12)/($10+$14))*100)}' | cut -f1-10,14- | awk -v idenVecGen=$idenVecGen '($17>=idenVecGen && $18>=idenVecGen)' > alignment.19.sam.A

#select columns and rearrange columns 
awk '{if ($2 ~ /chr/) print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $7 "\t" $6 "\t" $1 "\t "$8 "\t" $9 "\t" $10 "\t" $13 "\t" $15 "\t" $14 "\t" $16 "\t" $17 "\t" $18 "\t" $11 "\t" $12; else print $5 "\t" $7 "\t" $6 "\t" $2 "\t" $3 "\t" $4 "\t" $1 "\t "$8 "\t" $10 "\t" $9 "\t" $14 "\t" $16 "\t" $13 "\t" $15 "\t" $18 "\t" $17 "\t" $11 "\t" $12 }' alignment.19.sam.A  | awk '($11<=1 && $13<=1)' | sort -k8,8 -u > alignment.20.sam.A
###########################
echo " "
echo "Duplicate reads removal in process..."
cut -f7 alignment.20.sam.A | awk 'NR==FNR{tgts[$1]; next} $1 in tgts'  - alignment.2.sam.basic  > alignment.21.sam.A

$samtools view -bT $vecGenRef alignment.21.sam.A  > alignment.21.sam.A.bam
$samtools sort  alignment.21.sam.A.bam -o  alignment.21.sam.A.sorted.bam
$samtools rmdup alignment.21.sam.A.sorted.bam alignment.21.sam.A.sorted.nodup.bam
$samtools view  alignment.21.sam.A.sorted.nodup.bam  |cut -f1 |sort | uniq > alignment.21.sam.A.sorted.nodup.ids

awk 'NR==FNR{tgts[$1]; next} $7 in tgts' alignment.21.sam.A.sorted.nodup.ids alignment.20.sam.A > alignment.22.sam.A 

#split file into multiple based on vector
if [ -s "$outDir/alignment.22.sam.A" ]; then
rm -f *.split*
echo " "
echo "Clustering and annotation processing..."
sort -k4,4 -k1,1 -k2,2n alignment.22.sam.A  | awk '{print>$4".split"}'
sort -k12,12 *split | awk '{print>$4"."$12".split"}'

for f in *Hit.split ; do  sort -k1,1 -k2,2n $f | awk '{print $1"\t"$2"\t"$3"\t" $0}'  >> $f.txt  ; done

for f in *split.txt ; do  awk '{print $1"\t" $2"\t" $3}' $f | sort -k1,1 -k2,2 | uniq -c | awk '{print $2 "\t" $3"\t" $4"\t" $1}' | sort -k1,1 -k2,2 | awk -v clusterRange=5 -f $bin/vecGenClus.awk |tail -n+2  | sed 's/ /\t/g' >> $f.cluster  ; done

for f in *split.txt.cluster ; do awk '{print $1"\t"$2"\t"$3 "\t" $4 "\t" FILENAME}'  $f  | cut -d '.' -f1 | awk '{print $1"\t"$2"\t"$5 "\t" $4}'   >> $f.1  ; done
for f in *split.txt.cluster ; do awk '{print $1"\t"$2"\t"$4 "\t" FILENAME}'  $f  | cut -d '.' -f1 | awk '{print $1"\t"$2"\t"$4}'   >> $f.2  ; done

awk '{print $3 "\t" $5 "\t" $6 "\t" $2 "\t" $4 "\t" $7 "\t" $8 "\t" $9 "\t"$10 "\t" $11 "\t" $13 "\t"}' $annoTable | sed 's/ \+ /\t/g' | sed  '1d' > $outDir/anno1

bash $bin/vecGenAnn.sh $outDir $bedtools $annoTable $bin

awk '{print $1"@"$2"@"$4"@"$7"\t"$0}' alignment.22.sam.A > alignment.22.sam.A.ed
#for f in *split.txt.cluster.1.anno ; do  sort -k2,2 -u -t ','  $f >> $f.final0 ; done
for f in *split.txt.cluster.1.anno ; do awk -F ',' '{print $1"@"$2"@"$3"\t"$0}' $f  | sort -k1,1 -u | cut -f2-  >> $f.final0 ; done
for f in *split.txt.cluster.1.anno.final0 ; do awk -F ',' '{print $1"@"$2"@"$3}' $f >> $f.final1 ; done
for f in *split.txt.cluster.1.anno.final0 ; do awk -F ',' '{print $1"@"$2"@"$3"\t"$0}' $f | sort -k1,1 >> $f.final1.final3.final2 ; done

for f in *split ; do awk -F '\t' '{print $1"@"$2"@"$4"\t"$0}' $f >> $f.txt.cluster.1.anno.final0.final1.append ; done


for f in *split.txt.cluster.1.anno.final0.final1 ; do grep -Fwf $f $f.append | sort -k1,1 -u | sort -k1,1 >> $f.final3  ; done

for f in *split.txt.cluster.1.anno.final0.final1.final3 ; do paste $f $f.final2 >> $f.final3ed  ; done

echo " "
echo "Generating final result files for host-vector fusions..."

echo -e "Chr,GenomicPosition,StrandGenomic,VectorName,VectorPosition,StrandVector,ReadID,Sequence,SpanGenomic,SpanVector,Feature1_Genomic,Feature2_Vector,IdentityPercGenomic,IdentityPercVector,OverlapFusion,DistanceFusion,SequenceCount,RefSeqID,GeneStrand,GeneName,GeneLength,DisttoTSS,Upstream,Downstream,IntronExon" >  ISGenomeVector.csv

cat *final3ed | awk '{if($1==$20) print $0; else print "ERROR"}' | awk -F '\t' '{print $2 "\t" $3 "\t" $4"\t" $5 "\t" $6"\t" $7 "\t" $8"\t" $9 "\t" $10"\t" $11 "\t" $12"\t" $13 "\t" $14"\t" $15"\t"  $16"\t"  $17"\t" $18 "\t" $19 "\t" $21}' | sed 's/\t/,/g' | cut -d ',' -f1-18,22- | cut -d ',' -f1-10,12,14- | sort -t ',' -k19,19nr >> ISGenomeVector.csv
head -1 ISGenomeVector.csv > ISGenomeVector.UniqueGenome.csv
awk -F ',' '($11=="UniqueHit")' ISGenomeVector.csv >> ISGenomeVector.UniqueGenome.csv

head -1 ISGenomeVector.csv >  ISGenomeVector.NonUniqueGenome.csv
awk -F ',' '($11=="Non-uniqueHit")' ISGenomeVector.csv >> ISGenomeVector.NonUniqueGenome.csv


echo -e "Chr,GenomicPosition,StrandGenomic,VectorName,VectorPosition,StrandVector,ReadID,Sequence,SpanGenomic,SpanVector,OverlapFusion,DistanceFusion" >  ISGenomeVector.Unclustered.csv
less  alignment.22.sam.A | cut -f1-10,17- | sed 's/\t/,/g' >>  ISGenomeVector.Unclustered.csv
else
echo " "
echo "NO GENOME-VECTOR FUSIONS OR IS FOUND."
fi
#mv alignment.3.sam alignment_vecGenFusion.sam
echo " "
echo "     Host-Vector Fusion Analysis Finished"
echo "..........................................."
rm -f  ContAnalysis_WithoutContaminantsReadPairs.R1.fastq.gz ContAnalysis_WithoutContaminantsReadPairs.R2.fastq.gz  alignment_vecGenFusion.sam  ext* eext* vecToRemoveIDS anno* alignment.* *UniqueHit.split.txt* final* *.bed temp2.fastq temp1.fastq  file* *split  *txt* comb
echo "Host-Vector Fusion Analysis finished" >> log
date >> log
###########################
