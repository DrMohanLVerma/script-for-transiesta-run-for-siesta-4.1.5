#!/bin/bash
# =========================================================== #
# 	**** Script for scat and tbtrans run ( Positive-bias) #
#                     (for siesta-4.1.5)                      # 
#============================================================ #
# Please follow the steps : 				        #
# 1) modify this script file as per requirement of your       #
# (don't change the name of file as  well as system name).    #
# It is expected that all the binaries exe are linked with    #
# /usr/local/bin   					        #
# In presence of proper *.psf file for corresponding elements #
# run this script for transiesta and tbtrans run using command#
#	 $ sh mlv_posi_bias_scat_run                          #
# The calculation may take long time depending on size of     #
# system and number of nodes(in parallel run).                #
# The Value of current will be given in tbt.out     	        #
# After successfull run this script will give IV-curve for 	#
# Positive bias region						#
#=============================================================#
# Author: Mohan L Verma, Computational Nanomaterial           #  
# Research lab, Department of Applied Physics, FET,           #
# SSGI, Shri Shanakaracharya Technical Campus-Bhilai          # 
# (Chhattisgarh)  INDIA, www.drmlv.in                         #
#=============================================================
# Dont forget to give feedback : drmohanlv@gmail.com          #
#           April 03   ver:4.1.5    year: 2023                #
#=============================================================
> IV.dat
mkdir posi_bias
cd  posi_bias

cp ../../Elec/elec.TSHS .   # copy TSHS file from step-1

mkdir cont   # read the comment at the end of this script.

for i in `seq 0.0 0.1 1.2` 
do

cp -r cont scat$i
cd scat$i
cp ../../*.psf .
cp ../../Elec/elec.TSHS .

cp ~/SIESTA/siesta-v4.1.5/Util/TS/TBtrans/tbtrans .
 
cat > scat.fdf <<EOF

 SystemName    scat                                                                                                                     
SystemLabel   scat

# ---------------------------------------------------------------------------
# Lattice
# ---------------------------------------------------------------------------

LatticeConstant             1.00 Ang

%block LatticeVectors
     5.438063      0.000000      0.000000
     0.000000     20.324000      0.000000
     0.000000      0.000000     15.242090
%endblock LatticeVectors

# ---------------------------------------------------------------------------
# Species and Atoms
# ---------------------------------------------------------------------------

NumberOfSpecies        3
NumberOfAtoms         14

%block ChemicalSpeciesLabel
  1   1  H
  2   6  C
  3  79  Au
%endblock ChemicalSpeciesLabel

# ---------------------------------------------------------------------------
# Atomic Coordinates
# ---------------------------------------------------------------------------

AtomicCoordinatesFormat Ang

%block AtomicCoordinatesAndAtomicSpecies
 3.08373	10.03502	1.32255	3
3.08373	10.05698	3.96754	3
2.71903	10.32400	6.19985	2
4.89503	10.21800	6.36985	1
0.54303	10.21300	6.36985	1
3.94103	10.20600	6.88885	2
1.49703	10.20400	6.88885	2
1.49703	10.11800	8.30485	2
3.94103	10.11900	8.30485	2
0.54303	10.10700	8.82385	1
4.89503	10.11000	8.82485	1
2.71903	10.00000	8.99485	2
3.08373	10.03502	11.27455	3
3.08373	10.05698	13.91954	3
%endblock AtomicCoordinatesAndAtomicSpecies

# K-points

%block kgrid_Monkhorst_Pack
1   0   0   0.0
0   8   0   0.0
0   0   3  0.0
%endblock kgrid_Monkhorst_Pack

PAO.BasisSize    SZP
PAO.EnergyShift  0.005 Ry
==================================================
==================================================
# General variables

ElectronicTemperature  100 K 
MeshCutoff           350. Ry
xc.functional         LDA           # Exchange-correlation functional
xc.authors            CA 
SpinPolarized .false.
SolutionMethod Transiesta 

==================================================
==================================================
# SCF variables

DM.MixSCF1   T
MaxSCFIterations      300           # Maximum number of SCF iter
DM.MixingWeight       0.03          # New DM amount for next SCF cycle
DM.Tolerance          1.d-4         # Tolerance in maximum difference
DM.UseSaveDM          true          # to use continuation files
DM.NumberPulay         5
Diag.DivideAndConquer  no
Diag.ParallelOverK     yes

==================================================
==================================================
# MD variables

MD.FinalTimeStep 1
MD.TypeOfRun CG
MD.NumCGsteps     000
MD.UseSaveXV      .true.

==================================================
==================================================
# Output variables

WriteMullikenPop                1
WriteBands                      .false.
SaveRho                         .false.
SaveDeltaRho                    .false.
SaveHS                          .false.
SaveElectrostaticPotential      True 
SaveTotalPotential              no
WriteCoorXmol                   .true.
WriteMDXmol                     .true.
WriteMDhistory                  .false.
WriteEigenvalues                yes

==================================================
==================================================

TS.Voltage   $i eV
%block TS.ChemPots
  Left
  Right
%endblock TS.ChemPots

%block TS.ChemPot.Left
  mu V/2
  contour.eq
    begin
      c-Left
      t-Left
    end
%endblock TS.ChemPot.Left
%block TS.ChemPot.Right
  mu -V/2
  contour.eq
    begin
      c-Right
      t-Right
    end
%endblock TS.ChemPot.Right

TS.Elecs.Bulk true
TS.Elecs.DM.Update cross-terms
TS.Elecs.GF.ReUse true
%block TS.Elecs
  Left
  Right
%endblock TS.Elecs

%block TS.Elec.Left
  HS ./../../elec.TSHS
  chem-pot Left
  semi-inf-dir -a3
  elec-pos begin 1
  used-atoms 2
%endblock TS.Elec.Left

%block TS.Elec.Right
  HS ./../../elec.TSHS
  chem-pot Right
  semi-inf-dir +a3
  elec-pos end -1
  used-atoms 2
%endblock TS.Elec.Right

TS.Contours.Eq.Pole    2.50000 eV
%block TS.Contour.c-Left
  part circle
   from  -28.00000 eV + V/2 to -10. kT + V/2
    points 2
     method g-legendre
%endblock TS.Contour.c-Left
%block TS.Contour.t-Left
  part tail
   from prev to inf
    points 10
     method g-fermi
%endblock TS.Contour.t-Left
%block TS.Contour.c-Right
  part circle
   from  -28.00000 eV - V/2 to -10. kT - V/2
    points 2
     method g-legendre
%endblock TS.Contour.c-Right
%block TS.Contour.t-Right
  part tail
   from prev to inf
    points 10
     method g-fermi
%endblock TS.Contour.t-Right

TS.Elecs.Eta    0.0001000000 eV
%block TS.Contours.nEq
  neq
%endblock TS.Contours.nEq
%block TS.Contour.nEq.neq
  part line
   from -|V|/2 - 5 kT to |V|/2 + 5 kT
    delta 0.01 eV
     method mid-rule
%endblock TS.Contour.nEq.neq



# TBtrans options

TBT.T.Eig 3
TBT.Elecs.Eta    0.0000136058 eV

%block TBT.Contours
  neq
%endblock TBT.Contours

%block TBT.Contour.neq
  part line
   from   -3.00000 eV to    3.00000 eV
    delta    0.01200 eV
     method mid-rule
%endblock TBT.Contour.neq

# It is advised to define a device region of
# particular interest

EOF

mpirun -np 14 siesta.exe scat.fdf | tee  scat.out   # for scat run 


mpirun -np 14 tbtrans scat.fdf | tee  tbt.out    # for tbtrans run 


IV_Curve=`grep 'Left -> Right, V \[V] / I \[A]:' tbt.out | tail -1 | awk '{print $12}'`
echo $i '   '$IV_Curve >> ../IV.dat

cd ..
rm -rf cont 
mkdir cont

cp  ./scat$i/scat.TSDE ./scat$i/scat.TSHS ./scat$i/scat.DM  cont  # copy these files for continuation of the next bias step.



done

sed -i '1 i\#========================================================#' IV.dat
sed -i '2 i\#    This data is extracted using get_IV_curve_script    #' IV.dat
sed -i '3 i\#Author : Dr Mohan L Verma, Computational Nanomaterial   #' IV.dat
sed -i '4 i\#Research lab, Department of Applied Physics, FET,  SSGI #' IV.dat
sed -i '5 i\#Shri Shanakaracharya Technical Campus-Bhilai (CG) INDIA #' IV.dat
sed -i '6 i\#========================================================#' IV.dat
xmgrace IV.dat &

