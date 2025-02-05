#!/bin/bash

# Go to src, execute code, and move the respective files
make clean
make
cd ./bin
./create_e2c_purityntuple
echo "Purity Ntuple ready!"
./create_e2c_efficiencyntuple
echo "Efficiency Ntuple ready!"
# ./create_e2c_corrntuple
# echo "Corr. Data Ntuple ready!"
