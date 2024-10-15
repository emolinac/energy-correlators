#!/bin/bash

cd ./src

root -b -q macro_print_norme2c_data_jetpt.cpp
root -b -q macro_print_norme2c.cpp
root -b -q "macro_print_norme2c_data_jetpt.cpp(1)"
root -b -q "macro_print_norme2c.cpp(1)"