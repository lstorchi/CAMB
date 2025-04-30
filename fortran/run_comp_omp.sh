export filelist="*.dat"
if [ -d serial ]; then
  cd serial
  rm -f $filelist
  cd ../
else
  echo "serial directory does not exist"
  exit
fi
if [ -d omp ]; then
  cd omp
  rm -f $filelist
  cd ../
else
  echo "omp directory does not exist"
  exit
fi

if [ "$#" -ne 1 ]; then
  echo "No input file provided"
  exit
fi

if [ -f $1 ]; then
  cp ./serial/Makefile ./
  make clean ; make 
  time ./camb $1 
  cp -f $filelist ./serial

  cp ./omp/Makefile ./
  make clean ; make 
  time ./camb $1 
  cp -f $filelist ./omp/ 
else
  echo "Input file does not exist"
  exit
fi


