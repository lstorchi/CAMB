if [ -d GPU ]; then
  cd GPU
  rm -f *.dat
  cd ../
else
  echo "GPU directory does not exist"
  exit
fi
if [ -d CPU ]; then
  cd CPU
  rm -f *.dat
  cd ../
else
  echo "CPU directory does not exist"
  exit
fi

if [ "$#" -ne 1 ]; then
  echo "No input file provided"
  exit
fi

if [ -f $1 ]; then
 
  cp ./CPU/Makefile ./
  cd ../forutils ; make clean; cd - ; make clean ; make 
  ./camb $1 
  cp -f *.dat ./CPU/
  
  cp ./GPU/Makefile ./
  cd ../forutils ; make clean; cd - ; make clean ; make 
  ./camb $1 
  cp -f *.dat ./GPU/
 
  cd CPU/
  for name in * ; do echo $name; diff -b -B $name ../GPU/$name; done 
  cd ../
else
  echo "Input file does not exist"
  exit
fi


