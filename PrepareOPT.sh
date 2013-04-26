for file in *.pdb
do 
{
file1=${file%.*}
echo $file1
cp head.com $file1.com
cp $file tempt.txt
tail -n+3 tempt.txt > tempt1.txt
sed -i '$ d' tempt1.txt
#cut -d " " -f3,5-8 tempt1.txt > tempt2.txt
#awk -F, '{print $3,$5,$6,$7,$8}' tempt1.txt > tempt2.txt
#awk '{print $0"\t"$1"\t"$3"\t"$4"\t"$8"\t"$9"\t"$10"\t"substr($0, index($0,$10))}' tempt1.txt > tempt2.txt
awk -F" " '{print $3,$6,$7,$8}' tempt1.txt >> $file1.com
cat $file1.com | echo " " >> $file1.com
cat $file1.com | echo " " >> $file1.com
}
done
