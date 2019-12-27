outDir=$1
Ref=$2
aligner=$3

cd $outDir
echo " "
echo "Sensitive Rescue mode in process..."
echo " "

#Sensitive mode for no alternative tagged reads
awk '($6 ~ /S/)' alignmentSens.sam | grep -v 'XP' | awk '{print $0, gsub(/S/,"\t",$6) }' | awk '($15==1)' > ext1
grep 'XP:Z:' alignmentSens.sam | cut -f1 | uniq | grep -v -Fwf - ext1 > ext0
mv ext0 ext1 
awk '{ split($6,a , "S"); print $0 "\t" a[1]}' ext1 | awk '($16 ~ /M/)' |sort -k1,1 -u > ext1a
awk '{print $16}' ext1a |sed -e 's/[A-Z]/& /g' | awk '{print $NF}'  | paste ext1a - | awk -F " " '{print ">"$1 "\n" substr($10,(length($10)-$17+1),$17)}' > ext1b.fa
$aligner  mem -t 8 -T 15 -k 10 -w 10 $Ref ext1b.fa > ext1b.sam
#strand and XP
awk '($2=="0" || $2=="16" || $2="4")' ext1b.sam | grep -v '@SQ' | sort -k1,1 -u | cut -f1-14 |paste -d '\t' ext1a - | awk '($1==$17)' | awk '($18=="0" || $18=="16")' | awk '{if (($18=="16") && ($22 !~ /S/)) print $0 "\t" $19",""+"$20","$22"50S"",""NA"",""NA"";"; if (($18=="0") && ($22 !~ /S/)) print $0 "\t" $19",""-"$20",""50S"$22",""NA"",""NA"";"; if (($18=="16") && ($22 ~ /S/)) print $0 "\t" $19",""+"$20","$22",""NA"",""NA"";"; if (($18=="0") && ($22 ~ /S/)) print $0 "\t" $19",""-"$20","$22",""NA"",""NA"";" }' >  ext1c
awk '{print $0 "\t" $22}' ext1c | awk 'BEGIN{FS="\t"} {gsub("[^S]","",$31); print $0,length($31);}' OFS='\t'  | awk '($32=="0" || $33=="1")' | cut -f1-30 > ext1c1
mv ext1c1 ext1c
#now it is same as .3.sam file
awk '{print $1"\t" $2 "\t" $3 "\t" $4"\t" $5 "\t" $6 "\t" $7"\t" $8 "\t" $9 "\t" $10"\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" "XP:Z:"$31} ' ext1c  > ext1d
#add xp to the new aligned read - use this later
awk '{if ($2==147 || $2==83 || $2==145 || $2==81 || $2==113 || $2==177 ) print $0 "\t",  $2="-"; else print $0 "\t" "+" }' ext1c |  awk '{if ($22 !~ /S/ && $18=="16") print $17"\t" $18 "\t" $19 "\t" $20"\t" $21 "\t" $22"50S" "\t" $23"\t" $24 "\t" $25 "\t" $26"\t" $27 "\t" $28 "\t" $29 "\t" $30 "\t" "XP:Z:"$3","$32 $4","$6",""NA"",""NA"; if ($22 !~ /S/ && $18=="0") print $17"\t" $18 "\t" $19 "\t" $20"\t" $21 "\t" "50S"$22 "\t" $23"\t" $24 "\t" $25 "\t" $26"\t" $27 "\t" $28 "\t" $29 "\t" $30 "\t" "XP:Z:"$3","$32 $4","$6",""NA"",""NA"; if ($22 ~ /S/) print $17"\t" $18 "\t" $19 "\t" $20"\t" $21 "\t" $22 "\t" $23"\t" $24 "\t" $25 "\t" $26"\t" $27 "\t" $28 "\t" $29 "\t" $30 "\t" "XP:Z:"$3","$32 $4","$6",""NA"",""NA" } ' > ext1e
#SM
#awk '($6 ~ /S/)' alignmentSens.sam | grep -v 'XP' | awk '{print $0, gsub(/S/,"\t",$6) }' | awk '($15==1)' > eext1
awk '{ split($6,a , "S"); print $0 "\t" a[1]}' ext1 | awk '($16 !~ /M/)' | sort -k1,1 -u > eext1a
awk '{print $0  "\t" $16}' eext1a  | awk -F " " '{print ">"$1 "\n" substr($10,$1,$17)}' > eext1b.fa
$aligner  mem -t 8 -T 15 -k 10 -w 10 $Ref eext1b.fa > eext1b.sam
awk '($2=="0" || $2=="16" || $2="4")' eext1b.sam | grep -v '@SQ' | sort -k1,1 -u | cut -f1-14 |paste -d '\t' eext1a - | awk '($1==$17)' | awk '($18=="0" || $18=="16")' | awk '{if (($18=="16") && ($22 !~ /S/)) print $0 "\t" $19",""-"$20",""50S"$22",""NA"",""NA"";"; if (($18=="0") && ($22 !~ /S/)) print $0 "\t" $19",""+"$20","$22"50S"",""NA"",""NA"";"; if (($18=="16") && ($22 ~ /S/)) print $0 "\t" $19",""-"$20","$22",""NA"",""NA"";"; if (($18=="0") && ($22 ~ /S/)) print $0 "\t" $19",""+"$20","$22",""NA"",""NA"";" }'  >  eext1c
awk '{print $0 "\t" $22}' eext1c | awk 'BEGIN{FS="\t"} {gsub("[^S]","",$31); print $0,length($31);}' OFS='\t'  | awk '($32=="0" || $33=="1")' | cut -f1-30 > eext1c1
mv eext1c1 eext1c
#now it is same as .3.sam file
awk '{print $1"\t" $2 "\t" $3 "\t" $4"\t" $5 "\t" $6 "\t" $7"\t" $8 "\t" $9 "\t" $10"\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" "XP:Z:"$31} ' eext1c  > eext1d
#add xp to the new aligned read - use this later
awk '{if ($2==147 || $2==83 || $2==145 || $2==81 || $2==113 || $2==177 ) print $0 "\t",  $2="-"; else print $0 "\t" "+" }' eext1c |  awk '{if ($22 !~ /S/ && $18=="16") print $17"\t" $18 "\t" $19 "\t" $20"\t" $21 "\t""50S"$22 "\t" $23"\t" $24 "\t" $25 "\t" $26"\t" $27 "\t" $28 "\t" $29 "\t" $30 "\t" "XP:Z:"$3","$32 $4","$6",""NA"",""NA"; if ($22 !~ /S/ && $18=="0") print $17"\t" $18 "\t" $19 "\t" $20"\t" $21 "\t"$22"50S" "\t" $23"\t" $24 "\t" $25 "\t" $26"\t" $27 "\t" $28 "\t" $29 "\t" $30 "\t" "XP:Z:"$3","$32 $4","$6",""NA"",""NA"; if ($22 ~ /S/) print $17"\t" $18 "\t" $19 "\t" $20"\t" $21 "\t" $22 "\t" $23"\t" $24 "\t" $25 "\t" $26"\t" $27 "\t" $28 "\t" $29 "\t" $30 "\t" "XP:Z:"$3","$32 $4","$6",""NA"",""NA"} ' > eext1e

cat ext1d alignmentSens.sam eext1d > alignmentSens.sam.1
cat ext1e alignmentSens.sam.1 eext1e > alignmentSens.3.sam.ext

