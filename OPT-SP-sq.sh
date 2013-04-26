line="19"
sp="sp"
for file in *.out
do
{
echo $file
file1=${file%.*}
file2=$file1$sp
sed '/Input orientation/,/Distance/!d' $file|sed "s/tomic      Atomic             Coordinates (Angstroms)/#   /g"|grep -v 'Num'|grep -v 'Rot'|grep -v '\-\-\-'|grep -v 'S'|cut -c14-80|sed 's/	/       /g'|grep -v 'Distance'| grep -v Input > $file1.temp
tail --lines=$line $file1.temp > $file1.temp1
paste -d "" symbol.com $file1.temp1 > $file1.temp2
awk -F" " '{print $1,$4,$5,$6}' $file1.temp2 > $file1.temp3
cp head-sp.com $file2.com
cat $file1.temp3>>$file2.com
cat $file2.com | echo " " >> $file2.com
cat $file2.com | echo " " >> $file2.com
rm -rf $file1.temp
rm -rf $file1.temp1
rm -rf $file1.temp2
rm -rf $file1.temp3
}
done

