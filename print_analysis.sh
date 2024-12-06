#!/bin/bash

mother_dir=$(pwd)

cd ${mother_dir}/src-analysis/purity

root -b -q macro_print_pairpurity_rl.cpp
root -b -q macro_print_pairpurity_rl_jet_pt.cpp

cd ${mother_dir}/src-analysis/efficiency

root -b -q macro_print_pairefficiency_rl.cpp
root -b -q macro_print_pairefficiency_rl_jet_pt.cpp

cd ${mother_dir}/src-analysis/background

root -b -q macro_print_matching_fraction_e2c.cpp
root -b -q macro_print_matching_fraction.cpp
root -b -q macro_print_matching_fraction_jetpt.cpp

