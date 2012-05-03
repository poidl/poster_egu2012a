cd run53/
make
cd ../run54
make
cd ../run55
make
cd ../run56
make
cd ..
mpirun -n 2 ./run53/HIM < ./run53/temp_in
mpirun -n 2 ./run54/HIM < ./run54/temp_in
mpirun -n 2 ./run55/HIM < ./run55/temp_in
mpirun -n 2 ./run56/HIM < ./run56/temp_in
