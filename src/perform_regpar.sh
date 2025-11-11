#!/bin/bash

niter_syst_window=3
niter_unf=4
niter_jet_unf=4

# Check the following variable is the same on names.h
ext=".root"

cd ../src-analysis

for ((i=0 ; i<${niter_syst_window} ; i++)); do
        syst_niter=$((${niter_unf}+i-1))

        root -b -q "macro_print_histocorreec_rl_jetpt_weightpt.cpp(${syst_niter},${niter_jet_unf},\"--get-regpar\")"
done

