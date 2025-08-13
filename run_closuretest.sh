#!/bin/bash

# Go to src, execute code, and move the respective files
make clean
make
cd ./bin

echo "Creating corrections ntuple"
./create_eec_paircorrectionsntuple_ct
echo "Corrections ntuple ready!"

echo "Creating corr data ntuple"
./create_eec_corrntuple_paircorr_ct
echo "Corr. data ntuple ready!"


echo "Producing results"
cd ../src-analysis/
root -b -q "macro_print_fullcorreec_paircorr_2dunf_ct(4, true, true)"
