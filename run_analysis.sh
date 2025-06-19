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

echo "Creating pair purity ntuple"
./create_e2c_corrntuple_paircorr
echo "Pair purity ntuple ready!"

# echo "Creating jet purity ntuple"
# ./create_jet_purityntuple
# echo "Jet purity ntuple ready!"
# echo "Creating jet efficiency ntuple"
# ./create_jet_efficiencyntuple
# echo "Jet efficiency ntuple ready!"

echo "Creating corr data ntuple"
./create_e2c_corrntuple_paircorr
echo "Corr. data ntuple ready!"


# echo "Producing results"
# cd ../src-analysis/
# root -b -q macro_print_fullcorre2c_logbin.cpp
# root -b -q macro_print_fullcorre2c.cpp
