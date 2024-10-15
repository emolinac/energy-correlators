#!/bin/bash

cd ./src-corrections

root -b -q macro_print_pairpurity_rl.cpp
root -b -q macro_print_pairpurity_rl_jet_pt.cpp
root -b -q "macro_print_pairpurity_rl.cpp(1)"
root -b -q "macro_print_pairpurity_rl_jet_pt.cpp(1)"