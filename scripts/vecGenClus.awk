function cluster(chr,pos,str,num) {print chr,pos,str,num}
BEGIN{} {if((chr==$1)&&(($2-pos)*($2-pos)<=clusterRange*clusterRange)){tab[$1"_"$2]+=$4; if(tab[$1"_"$2]>max){max=tab[$1"_"$2];chr=$1; pos=$2; str=$3} num+=$4} else {cluster(chr,pos,str,num); tab[$1"_"$2]+=$4; chr=$1; pos=$2; str=$3; num=$4; max=$4}} END {if(num>=countThr) cluster(chr,pos,str,num)}
