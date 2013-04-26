for file in *.inp
do
{
sed -i '/number 1/a\  resnumbers 2' $file
}
done
