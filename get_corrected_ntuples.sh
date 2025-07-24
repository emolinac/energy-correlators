#!/bin/bash

# Go to src, execute code, and move the respective files
make clean
make
cd ./bin

# ./create_e2c_paircorrectionsntuple
# ./create_e2c_paircorrectionsntuple_ct
# ./create_e2c_paircorrectionsntuple_prior
# ./create_e2c_paircorrectionsntuple_jer
# ./create_e2c_paircorrectionsntuple_jes
# ./create_e2c_hadroncorrectionsntuple
# ./create_e2c_mc_ntuple

# ./create_jes_jer_ntuple
# ./create_jet_efficiencyntuple
# ./create_jet_purityntuple

./create_e2c_corrntuple
./create_e2c_corrntuple_paircorr
./create_e2c_corrntuple_paircorr_ct
./create_e2c_corrntuple_paircorr_prior
./create_e2c_corrntuple_paircorr_jer
./create_e2c_corrntuple_paircorr_jes

