#!/bin/bash
g++ -O3 -std=c++11 -I. -o prepare.bin prepare.cpp
./prepare.bin H.txt $(ls O*.txt)
g++ -O3 -std=c++11 -o PMRQMC.bin PMRQMC.cpp
./PMRQMC.bin
