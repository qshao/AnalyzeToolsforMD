for file in *.pdb
do 
{
file1=${file:0:7}
echo $file1
echo 1 | /usr/share/gromacs45/bin/pdb2gmx -f $file1.pdb -o $file1.gro -p $file1.top -water tip4p
}
done
