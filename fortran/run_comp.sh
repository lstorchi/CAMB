if [ -d GPU ]; then
  cd GPU
  rm -f test_params.ini s_points.txt Source_q.txt s_npoints.txt iv_sourcessteps.txt ddSource_q.txt test_scalCls.dat test_scalarCovCls.dat test_lensedCls.dat test_lenspotentialCls.dat 
  cd ../
else
  echo "GPU directory does not exist"
  exit
fi
if [ -d CPU ]; then
  cd CPU
  rm -f test_params.ini s_points.txt Source_q.txt s_npoints.txt iv_sourcessteps.txt ddSource_q.txt test_scalCls.dat test_scalarCovCls.dat test_lensedCls.dat test_lenspotentialCls.dat 
  cd ../
else
  echo "CPU directory does not exist"
  exit
fi

if [$# -ne 1]; then
  echo "No input file provided"
  exit
fi

if [ -f $1 ]; then
  cp ./GPU/Makefile ./
  cd ../forutils ; make clean; cd - ; make clean ; make 
  ./camb $1 
  mv test_params.ini s_points.txt Source_q.txt s_npoints.txt iv_sourcessteps.txt ddSource_q.txt test_scalCls.dat test_scalarCovCls.dat test_lensedCls.dat test_lenspotentialCls.dat ./GPU/
  
  cp ./CPU/Makefile ./
  cd ../forutils ; make clean; cd - ; make clean ; make 
  ./camb $1 
  mv test_params.ini s_points.txt Source_q.txt s_npoints.txt iv_sourcessteps.txt ddSource_q.txt test_scalCls.dat test_scalarCovCls.dat test_lensedCls.dat test_lenspotentialCls.dat ./CPU/
  
  cd CPU/
  for name in * ; do echo $name; diff -b -B $name ../GPU/$name; done 
  cd ../
else
  echo "Input file does not exist"
  exit
fi


