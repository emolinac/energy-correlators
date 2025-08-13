#!/bin/bash

# Go to src, execute code, and move the respective files
make clean
make
cd ./bin

echo "Creating purity ntuple"
./create_eec_purityntuple
echo "Purity ntuple ready!"
echo "Creating efficiency ntuple"
./create_eec_efficiencyntuple
echo "Efficiency ntuple ready!"

echo "Creating pair purity ntuple"
./create_eec_pairpurityntuple
echo "Pair purity ntuple ready!"
echo "Creating pair efficiency ntuple"
./create_eec_pairefficiencyntuple
echo "Pair efficiency ntuple ready!"

# echo "Creating jet purity ntuple"
# ./create_jet_purityntuple
# echo "Jet purity ntuple ready!"
# echo "Creating jet efficiency ntuple"
# ./create_jet_efficiencyntuple
# echo "Jet efficiency ntuple ready!"

echo "Creating corr data ntuple"
./create_eec_corrntuple_paircorr
echo "Corr. data ntuple ready!"


# echo "Producing results"
# cd ../src-analysis/
# root -b -q macro_print_fullcorreec.cpp
# root -b -q macro_print_fullcorreec.cpp
