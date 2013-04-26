sol="-sol"
app="1"
for file in *.gro
do 
{
file1=${file:0:7}
file2=$file1$sol
echo $file2
/usr/share/gromacs45/bin/genbox -cs tip4p.gro -cp $file1.gro -box 2.5 2.5 2.5  -o $file2$app.gro
echo q | /usr/share/gromacs45/bin/make_ndx -f $file2$app.gro -o $file2.ndx
}
done
