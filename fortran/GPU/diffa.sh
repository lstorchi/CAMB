for name in *.dat
do 
  echo $name
  diff $name ../GPU/$name 
done
