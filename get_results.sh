#!/bin/bash

cd ./src-analysis

root -b -q macro_print_fullcorreec.cpp
root -b -q macro_print_fullcorreec_paircorr_2dunf.cpp
root -b -q macro_print_fullcorreec_paircorr_2dunf_ct_niter.cpp
root -b -q macro_print_fullcorreec_paircorr_2dunf_jer.cpp
root -b -q macro_print_fullcorreec_paircorr_2dunf_jes.cpp
root -b -q macro_print_fullcorreec_paircorr_2dunf_prior.cpp
root -b -q macro_print_fullcorreec_paircorr_2dunf_shapect_niter.cpp

root -b -q macro_print_deviation_from_nominal_ct.cpp
root -b -q macro_print_deviation_from_nominal_shapect.cpp
root -b -q "macro_print_deviation_from_nominal.cpp(true,true,0)"
root -b -q "macro_print_deviation_from_nominal.cpp(true,true,2)"
root -b -q "macro_print_deviation_from_nominal.cpp(true,true,3)"
root -b -q "macro_print_deviation_from_nominal.cpp(true,true,4)"

root -b -q macro_print_fullcorreec_paircorr_2dunf_incsyst.cpp

cd ..