#!/bin/bash
make clean
make all -j
out_folder=/scratch1/07887/zipeng/Magnus_runs/
cp Main_ScalarField3d.Linux.64.mpiicpc.ifort.OPTHIGH.MPI.OPENMPCC.Skylake.Intel2018.ex $out_folder