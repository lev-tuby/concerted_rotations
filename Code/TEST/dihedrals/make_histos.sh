#!/bin/bash
S=100
awk -v s=$S 'NR%s==0 {print}' dihedral*.dat | ~/Scripts/histo_1d.awk col=13 min=-3.14157 max=3.14157 nbin=25  > phi_13.dat
awk -v s=$S 'NR%s==0 {print}' dihedral*.dat | ~/Scripts/histo_1d.awk col=14 min=-3.14157 max=3.14157 nbin=25  > psi_14.dat
awk -v s=$S 'NR%s==0 {print}' dihedral*.dat | ~/Scripts/histo_1d.awk col=15 min=-3.14157 max=3.14157 nbin=25  > phi_15.dat
awk -v s=$S 'NR%s==0 {print}' dihedral*.dat | ~/Scripts/histo_1d.awk col=16 min=-3.14157 max=3.14157 nbin=25  > psi_16.dat
awk -v s=$S 'NR%s==0 {print}' dihedral*.dat | ~/Scripts/histo_1d.awk col=17 min=-3.14157 max=3.14157 nbin=25  > phi_17.dat
awk -v s=$S 'NR%s==0 {print}' dihedral*.dat | ~/Scripts/histo_1d.awk col=18 min=-3.14157 max=3.14157 nbin=25  > psi_18.dat
awk -v s=$S 'NR%s==0 {print}' dihedral*.dat | ~/Scripts/histo_1d.awk col=19 min=-3.14157 max=3.14157 nbin=25  > phi_19.dat
awk -v s=$S 'NR%s==0 {print}' dihedral*.dat | ~/Scripts/histo_1d.awk col=20 min=-3.14157 max=3.14157 nbin=25  > psi_20.dat

#~/Scripts/histo_1d.awk col=3 min=-3.14157 max=3.14157 nbin=50 dihedral*.dat > phi_3.dat
#~/Scripts/histo_1d.awk col=4 min=-3.14157 max=3.14157 nbin=50 dihedral*.dat > psi_4.dat
#~/Scripts/histo_1d.awk col=5 min=-3.14157 max=3.14157 nbin=50 dihedral*.dat > phi_5.dat
#~/Scripts/histo_1d.awk col=6 min=-3.14157 max=3.14157 nbin=50 dihedral*.dat > psi_6.dat
#~/Scripts/histo_1d.awk col=7 min=-3.14157 max=3.14157 nbin=50 dihedral*.dat > phi_7.dat
#~/Scripts/histo_1d.awk col=8 min=-3.14157 max=3.14157 nbin=50 dihedral*.dat > psi_8.dat
