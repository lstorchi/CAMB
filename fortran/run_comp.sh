cp ./GPU/Makefile ./
cd ../forutils ; make clean; cd - ; make clean ; make -j 8
./camb  small.ini 
mv test_params.ini s_points.txt Source_q.txt s_npoints.txt iv_sourcessteps.txt ddSource_q.txt test_scalCls.dat test_scalarCovCls.dat test_lensedCls.dat test_lenspotentialCls.dat ./GPU/

cp ./CPU/Makefile ./
cd ../forutils ; make clean; cd - ; make clean ; make -j 8
mv test_params.ini s_points.txt Source_q.txt s_npoints.txt iv_sourcessteps.txt ddSource_q.txt test_scalCls.dat test_scalarCovCls.dat test_lensedCls.dat test_lenspotentialCls.dat ./CPU/

cd CPU/
for name in * ; do echo $name; diff -b -B $name ../GPU/$name; done 
cd ../
