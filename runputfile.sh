for file in *.top
do 
{
file1=${file:0:3}
file2=${file:0:7}
mkdir $file1
mkdir $file1/$file2
cp $file2*.* $file1/$file2
}
done
