#!/bin/bash

cd ./src-corrections

root -b -q macro_print_pairefficiency_rl.cpp
root -b -q macro_print_pairefficiency_rl_jet_pt.cpp
root -b -q "macro_print_pairefficiency_rl.cpp(1)"
root -b -q "macro_print_pairefficiency_rl_jet_pt.cpp(1)"