#!/bin/bash

# Go to src, execute code, and move the respective files
make clean
make
cd ./bin

# ./create_eec_paircorrectionsntuple
# ./create_eec_paircorrectionsntuple_ct
# ./create_eec_paircorrectionsntuple_prior
# ./create_eec_paircorrectionsntuple_jer
# ./create_eec_paircorrectionsntuple_jes
# ./create_eec_hadroncorrectionsntuple
# ./create_eec_mc_ntuple

# ./create_jes_jer_ntuple
# ./create_jet_efficiencyntuple
# ./create_jet_purityntuple

./create_eec_corrntuple
./create_eec_corrntuple_paircorr
./create_eec_corrntuple_paircorr_ct
./create_eec_corrntuple_paircorr_prior
./create_eec_corrntuple_paircorr_jer
./create_eec_corrntuple_paircorr_jes

