outDir=$1
file1=$2
file2=$3
qual=$4
lenPer=$5
adapter1=$6
adapter2=$7
trimmer=$8

cd $outDir

echo "Analysis started" > log
date >> log
echo " "
echo "Quality Control in process"
echo " "
echo "$trimmer -x $adapter1 -y $adapter2 -q $qual -l $lenPer -o $outDir/file $file1 $file2  "
echo " "
$trimmer -x  $adapter1 -y  $adapter2 -q $qual -l $lenPer -o $outDir/file $file1 $file2
echo " "
echo "Quality Control is finished"
echo " "
