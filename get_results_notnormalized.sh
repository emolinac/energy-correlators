#!/bin/bash

cd ./src-analysis

root -b -q -l macro_print_fullcorreec.cpp
root -b -q -l macro_print_fullcorreec_paircorr_2dunf.cpp
root -b -q -l "macro_print_fullcorreec_paircorr_2dunf_ct_niter(4, 10, true, true, false)"
root -b -q -l macro_print_fullcorreec_paircorr_2dunf_jer.cpp
root -b -q -l macro_print_fullcorreec_paircorr_2dunf_jes.cpp
root -b -q -l macro_print_fullcorreec_paircorr_2dunf_prior.cpp
root -b -q -l "macro_print_fullcorreec_paircorr_2dunf_shapect_niter(4, 1, true, true, false)"

root -b -q -l macro_print_deviation_from_nominal_ct.cpp
root -b -q -l macro_print_deviation_from_nominal_shapect.cpp
root -b -q -l "macro_print_deviation_from_nominal.cpp(false,true,0)"
root -b -q -l "macro_print_deviation_from_nominal.cpp(false,true,2)"
root -b -q -l "macro_print_deviation_from_nominal.cpp(false,true,3)"
root -b -q -l "macro_print_deviation_from_nominal.cpp(false,true,4)"

root -b -q -l macro_print_fullcorreec_paircorr_2dunf_incsyst.cpp

cd ..