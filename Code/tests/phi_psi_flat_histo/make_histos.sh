#!/bin/bash

## set how many steps
NUM=10000000
## dihedrals are printed every 50 iterations ... so each line corespond to 50 iterations
LINES=$(bc <<< "${NUM} / 50")

head -n ${LINES} dihedrals.dat > tmp.dat

PLOT="plot [][0:]"
PLOT_PHI="plot [][0:]"
PLOT_PSI="plot [][0:]"

for resi in $(seq 0 1 17); do
    ./histo_1d.awk col=$((resi * 2 + 1)) min=-3.14157 max=3.14157 nbin=50 tmp.dat > phi_${resi}.dat
    PLOT_PHI=$PLOT_PHI" 'phi_${resi}.dat' u 1:2 w boxes t '', "

    ./histo_1d.awk col=$((resi * 2 + 2)) min=-3.14157 max=3.14157 nbin=50 tmp.dat > psi_${resi}.dat
    PLOT_PSI=$PLOT_PSI" 'psi_${resi}.dat' u 1:2 w boxes t '', "

    PLOT=$PLOT" 'phi_${resi}.dat' u 1:2 w boxes t '', 'psi_${resi}.dat' u 1:2 w boxes t '', "
done

PLOT=${PLOT%%??}
PLOT_PHI=${PLOT_PHI%%??}
PLOT_PSI=${PLOT_PSI%%??}

gnuplot << EOF
set term pdf enhanced
set output 'phi_dihedrals.pdf'
$PLOT_PHI
EOF

gnuplot << EOF
set term pdf enhanced
set output 'psi_dihedrals.pdf'
$PLOT_PSI
EOF

gnuplot << EOF
set term pdf enhanced
set output 'dihedrals.pdf'
$PLOT
EOF

