outDir=${1}
bedtools=${2}
cd $outDir
for f in *split.txt.cluster.1
do
	sed 's/^ *//g'  $f | awk '{print ($1"\t"$2"\t"$2+1"\t"$3"\t"$4)}' > file2
	$bedtools intersect -a anno1 -b file2 -wa -wb | awk '{print $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11,$17 }' > file3a 
	awk '{ n14 = split($14, t14, ",");n15 = split($15, t15, ","); for (i = 0; ++i < n14;) {print $1,t14[i],t15[i], $10=="+" ? "Exon" i : "Exon" n14-i, $10} }' file3a | awk '{print $1"\t"$2"\t"$3"\t"$4"\t""0""\t"$5}' > fileexons.bed
	$bedtools intersect -a fileexons.bed  -b file3a  -wb > file3b
	awk '{n14 = split($14, t14, ",");n15 = split($15, t15, ",");for (i = 0; ++i < n14-1;) {m14= n14-1; print $1,t15[i],t14[i+1], $10=="+" ? "Intron" i : "Intron"  m14-i, $10}}' file3a | awk '{print $1"\t"$2"\t"$3"\t"$4"\t""0""\t"$5}' > fileintrons.bed
        $bedtools intersect -a fileintrons.bed -b file3a -wb > file3c
        cat file3c file3b | sort -k1,1 -k2,2n -u > file3d
	awk '{print $1 "\t" $2 "\t" $3 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16 "\t" $17 "\t" $18 "\t" $19 "\t" $20 "\t" $21 "\t" $22 "\t" $23 "\t" $4}' file3d > file3e
	cut -f1,2 file3e > file3eids
	cut -f1,2 file2 > file2ids
	grep -v -Fwf file3eids file2ids > fileNotIntExtids
	grep -Fwf fileNotIntExtids file2  > fileNotIntExt
	$bedtools closest -a fileNotIntExt -b anno1  -t  first | sed 's/ /\t/g'  > fileNotIntExt1
	cat fileNotIntExt1 file3e | awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5"@"$17 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $13 "\t" $14 "\t" $15 "\t" $16}' > comb
	#RefSeqID,GeneStrand,GeneName,GeneLength,DistTSS,Upstream,Downstream,IntronExon
	awk '{print ($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$8-$7) }' comb > file4
	awk '{ if ($2>=$7 && $2<=$8 && $10 ~ /+/) print ($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$2-$7-1) ; if ($2>=$7 && $2<=$8 && $10 ~ /-/) print ($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$8-$2); if ($2<$7 || $2>$8) print ($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$10)}' file4  > file5
	awk '{if($2<$7 && $10 ~ /+/) print ($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$7-$2+1); if ($2>$8 && $10 ~ /-/) print ($1"\t"$2"\t" $3"\t" $4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$2-$8) ; if ($2>=$7 && $2<=$8 && $10 ~ /+/ || $2>=$7 && $2<=$8 && $10 ~ /-/) print ($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$10); if ($2>$8 && $10 ~ /+/) print ($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18,$10); if ($2<$7 && $10~ /-/) print ($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$10)}' file5  > file6
	awk '{ if ($2>$8 && $10   ~ /+/) print ($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$2-$8); if ($2<$7 && $10 ~ /-/) print ($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$7-$2+1); if ($2>=$7 && $2<=$8 && $10 ~ /+/ || $2>=$7 && $2<=$8 && $10 ~ /-/) print ($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$10) ; if ($2<$7 && $10  ~ /+/) print ($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$10) ; if ($2>$8 && $10  ~ /-/) print ($1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$10)  } ' file6 > file7
	awk '{print $1 "\t" $2 "\t" $4 "\t" "count" "\t" $9 "\t" $10 "\t" $16 "\t" $17 "\t" $18 "\t" $19 "\t" $20 "\t" $5}' file7 | sed 's/@/\t/g' > file8
	awk '(($10=="-" && $11=="-") || ($10=="+" && $11=="+"))' file8 | awk '{print $1","$2","$3","$12","$5","$6","$7","$8","$9","" "","" "","$13}' > file9
	awk '(($9=="-" || $9=="+" ))' file8 |awk '($10=="-" || $10=="+")'|  awk '{print $1","$2","$3","$12","$5","$6","$7","$8","" "","" "","$11","""}' | sed 's/$@//g' >> file9
	awk '(($9=="-" || $9=="+" ))' file8 |awk '($11=="-" || $11=="+")'|  awk '{print $1","$2","$3","$12","$5","$6","$7","$8","" "","$10","" "","""}'  >> file9
cut -d ',' -f1,2 file9 | sed 's/,/\t/g' | grep -v -Fwf - file2 | awk '{print $1, $2, $4, $5, "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA"}' OFS="," | cat file9 - > $f.anno
done

#This sript contains modified version of a part of script from published method Afzal S et al. Mol Ther Nucleic Acid. 2017;6(March):133-139. doi:10.1016/j.omtn.2016.12.001
