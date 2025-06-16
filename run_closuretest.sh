#!/bin/bash

# Go to src, execute code, and move the respective files
make clean
make
cd ./bin

echo "Creating pair purity ntuple"
./create_e2c_pairpurityntuple_ct
echo "Pair purity ntuple ready!"
echo "Creating pair efficiency ntuple"
./create_e2c_pairefficiencyntuple_ct
echo "Pair efficiency ntuple ready!"

echo "Creating corr data ntuple"
./create_e2c_corrntuple_paircorr_ct
echo "Corr. data ntuple ready!"


echo "Producing results"
cd ../src-analysis/
root -b -q "macro_print_fullcorre2c_paircorr_2dunf_ct(4, true, true)"
