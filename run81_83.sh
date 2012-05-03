cd run81/
make
cd ../run82
make
cd ../run83
make
cd ..
mpirun -n 2 ./run81/HIM < ./run81/temp_in
mpirun -n 2 ./run82/HIM < ./run82/temp_in
mpirun -n 2 ./run83/HIM < ./run83/temp_in

