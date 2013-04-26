sol="-sol"
a1="1"
a2="2"
for file in *.top
do 
{
file1=${file:0:3}
file2=${file:0:7}
file3=$file2$sol
cp mdp1.mdp $file1/$file2/$file3$a1.mdp
cp mdp2.mdp $file1/$file2/$file3$a2.mdp
cp runAs.sh $file1/$file2/runAs.sh
sed -i -e "s/file/$file3/g" $file1/$file2/runAs.sh
}
done
