export filelist="out ddScaledSrc.txt test_params.ini iv_sourcessteps.txt ScaledSrc.txt s_points.txt Source_q.txt s_npoints.txt ddSource_q.txt test_scalCls.dat test_scalarCovCls.dat test_lensedCls.dat test_lenspotentialCls.dat"
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
  cp ./omp/Makefile ./
  make clean ; make 
  time ./camb $1 > out
  mv $filelist ./omp/ 
  
  cp ./serial/Makefile ./
  make clean ; make 
  time ./camb $1 > out
  mv $filelist ./serial
else
  echo "Input file does not exist"
  exit
fi


