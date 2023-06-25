#!/bin/bash
#
# This program is introduced in the paper:
# Lev Barash, Arman Babakhani, Itay Hen, A quantum Monte Carlo algorithm for arbitrary spin-1/2 Hamiltonians (2023).
#
# This program is licensed under a Creative Commons Attribution 4.0 International License:
# http://creativecommons.org/licenses/by/4.0/
#
#
g++ -O3 -std=c++11 -o prepare.bin prepare.cpp
./prepare.bin H.txt $(ls O*.txt  2> /dev/null)
g++ -O3 -std=c++11 -o PMRQMC.bin PMRQMC.cpp
./PMRQMC.bin
