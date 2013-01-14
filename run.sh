#!/bin/bash
./quantfdtd < freepar > plot_gnu 
gnuplot free_part.gnu 
gnuplot fft.gnu
