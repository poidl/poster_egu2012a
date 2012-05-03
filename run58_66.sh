cd run58/
make
cd ../run60
make
cd ../run61
make
cd ../run62
make
cd ../run63
make
cd ../run64
make
cd ../run66
make
cd ..
mpirun -n 2 ./run58/HIM < ./run58/temp_in
mpirun -n 2 ./run60/HIM < ./run60/temp_in
mpirun -n 2 ./run61/HIM < ./run61/temp_in
mpirun -n 2 ./run62/HIM < ./run62/temp_in
mpirun -n 2 ./run63/HIM < ./run63/temp_in
mpirun -n 2 ./run64/HIM < ./run64/temp_in
mpirun -n 2 ./run66/HIM < ./run66/temp_in

