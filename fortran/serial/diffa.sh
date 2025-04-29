for name in * 
do 
  echo $name 
  diff $name ../omp/$name 
done
