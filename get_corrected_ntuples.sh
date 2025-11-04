#!/bin/bash

# Go to src, execute code, and move the respective files
make clean
make
cd ./bin

# ./create_ntuple_reco2truth_match
# ./create_ntuple_reco2truth_match_ct
# ./create_ntuple_reco2truth_match_prior
# ./create_ntuple_reco2truth_match_jer
# ./create_ntuple_reco2truth_match_jes
# ./create_eec_hadroncorrectionsntuple
# ./create_eec_mc_ntuple

# ./create_jes_jer_ntuple
# ./create_jet_ntuple_truth2reco_match
# ./create_jet_purityntuple

./create_eec_corrntuple
./create_eec_corrntuple_paircorr
./create_eec_corrntuple_paircorr_ct
./create_eec_corrntuple_paircorr_prior
./create_eec_corrntuple_paircorr_jer
./create_eec_corrntuple_paircorr_jes

