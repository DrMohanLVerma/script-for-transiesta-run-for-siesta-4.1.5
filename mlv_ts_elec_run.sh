#!/bin/bash
# ===========================================================#
# 	**** Script for electrode run ****        	       # 
#============================================================#
# Please follow the steps : 				       #
# 1)make the electrode input file elec.fdf  for your system  #
# (don't change the name of file as  well as system name).   #
# It is expected that all the binaries exe are linked with   #
# /usr/local/bin   					       #
# In presence of proper *.psf file for corresponding elements#
# run this script for lectrode using command                 #  
#	 $ mlv_ts_elec_run                                   #
# The calculation should complete in a few minutes and will  #
#  generate a elec.TSHS file.		        	       #
#=============================================================
#=============================================================#
# Author: Mohan L Verma, Computational Nanomaterial           #  
# Research lab, Department of Applied Physics, FET,           #
# SSGI, Shri Shanakaracharya Technical Campus-Bhilai          # 
# (Chhattisgarh)  INDIA, www.drmlv.in                         #
#=============================================================
# Dont forget to give feedback : drmohanlv@gmail.com          #
#           April 03   ver:4.1.5    year: 2023                #
#=============================================================
 
echo "Electrode Calculation"

mkdir Elec

cd Elec

cp ../elec.fdf .
cp ../*.psf .
 

mpirun -np 14 siesta.exe --electrode elec.fdf | tee elec.out 


 
cp elec.out  ../ts_elec.out
#
# Go back to base directory


#
cd ..

 