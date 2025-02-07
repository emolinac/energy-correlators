#!/bin/bash

# Go to src, execute code, and move the respective files
make clean
make
cd ./bin
echo "Creating purity ntuple"
./create_e2c_purityntuple
echo "Purity ntuple ready!"
echo "Creating efficiency ntuple"
./create_e2c_efficiencyntuple
echo "Efficiency ntuple ready!"
echo "Creating corr data ntuple"
./create_e2c_corrntuple
echo "Corr. data ntuple ready!"
