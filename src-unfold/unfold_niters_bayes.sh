#!/bin/bash

for ((num = 1; num <= 5; num++))
do
    # root -b -q "macro_print_unfold_bayes_eta.cpp(${num},30,50)"
    root -b -q "macro_print_unfold_bayes_rl.cpp(${num},20,30)"
done