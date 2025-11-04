#!/bin/bash

niter_jer=10
niter_unf=4
niter_jet_unf=4

# Check the following variable is the same on names.h
std_output_name="histos_3dpaircorr_rl_jetpt_weightpt_eec_jer.root"
name="histos_3dpaircorr_rl_jetpt_weightpt_eec_jer"
ext=".root"

for ((i=0 ; i<${niter_jer} ; i++)); do
        cd ../bin

        ./create_ntuple_reco2truth_match --get-jer

        ./create_correec_histo3dpaircorr_rl_jetpt_weightpt --get-jer

        cd ../output-files
        mv ${std_output_name} ${name}_${i}${ext}
done

cd ../src-analysis

for ((i=0 ; i<${niter_jer} ; i++)); do
        root -b -q "macro_print_histocorreec_rl_jetpt_weightpt.cpp(${niter_unf},${niter_jet_unf},\"--get-jer\",${i})"

        outfile_name="histos_eec_3dcorr_rl_jetpt_weightpt_niter${niter_unf}_niterjet${niter_jet_unf}--get-jer"

        cd ../output-files
        mv ${outfile_name}${ext} ${outfile_name}_${i}${ext}

        cd ../src-analysis
done

