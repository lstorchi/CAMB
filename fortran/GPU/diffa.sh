for name in * 
do 
  echo $name
  diff $name ../GPU/$name 
done
