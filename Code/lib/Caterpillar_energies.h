#ifndef __CAT_ENER__
#define __CAT_ENER__
/**
 * @file
 * @brief File contain function for energy calculation
 */

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include "Caterpillar.h"

/** @brief Minimal distance between Ca atoms (hard-core repulsion cut-off). */
#define CATENR_SAW_radius   4.0000000

/** @brief Distance at which interaction energy is half of the well depth @todo add picture of square well like potential. */
#define CATENR_CaCa_RANGE   12.0000000

/** @brief Longest distance at which Ca atoms are assumed to interact. */
#define CATENR_CaCa_CUTOFF  17.0000000

/** @brief Cut-off distance for distance between H and O in hydrogen bond. */
#define CATENR_HB_CUTOFF		6.0000000

/** @brief Scaling factor for hydrogen bond interactions. 3.1*8*0.56, with truncation er. */
#define CATENR_HB_PREF   	  13.888

/** @brief Scaling factor for Ca Ca attractive interactions. 1-0.56, with truncation er. */
#define CATENR_CaCa_PREF    0.4400000000000000

/** @brief Not sure but assume default Hydropathy scale from aapot. */
#define CATENR_Hydropathy   0.014655

/** @brief Hydropathy scaling factor. CaCa_PREF*Hydropathy */
#define CATENR_WAT_PREF     0.00644819999999999926

/** @brief Value substituted for infinite interaction due to hard-core overlap. */
#define CATENR_INFTY 100000000000
/** @brief Definition of zero energy. */
#define CATENR_ZERO_EN 0

/**
 * @brief Function calculate hard-core repulstion energy between two atoms
 */
double CATENR_Saw 	    (double *r_atom_i, double *r_atom_j, double saw);

/**
 * @brief Function calculate hard-core repulstion energy from distance and cut-off
 */
double CATENR_Saw_eq 	(double r, double saw);

/**
 * @brief Function calculate attractive part of Catepillar interaction
 *
 * Function calculate continuous square well like potential.
 * By multipliing this value with Ca-Ca interaction energy we get whole Ca-Ca interaction.
 */
double CATENR_bond 			(double r);

/**
 * @brief Function calculate attractive part of Catterpillar interaction
 *
 * Function take atom types of both Ca atoms, interaction matrix M and preffactor.
 * bond is calculated via #CATENR_bond.
 */
double CATENR_Ca_Ca_eq	(int type_i, int type_j, double bond, double **M, double pref);

/**
 * @brief All in one function for calculation of attractive part of Caterpillar model
 *
 * LT tested OK 11/07/17 
 * @note Is it right that prefactor pref is being used even if LINK_DIAGNOSE macro being invoked?
 * @todo remove n variable ... not used at all.
 */
double CATENR_Ca_Ca 		(double *r_Ca_i, int type_i, double *r_Ca_j, int type_j, int n, double **M,double pref);

/**
 * @brief Calculate H-bond interaction from positions of N H and O atoms
 *
 * LT Tested OK (?) 12/07/17. There is still a small discrepancy that I suppose is caused by pdb errors. must check with the binary file.
 * @todo Since prefactor is given as macro it migh be removed from parameters ... unless same function being used for optimization purposes...
 */
double CATENR_HB  			(double *restrict r_H_i, double *restrict r_N_i, double *restrict r_O_j, double *restrict r_C_j, double pref);

/**
 * @brief Calculate H-bond interaction from positions of N H and O atoms
 *
 * LT Tested OK (?) 12/07/17. There is still a small discrepancy that I suppose is caused by pdb errors. must check with the binary file.
 * @note This function in comparison to #CATENR_HB do not check if atoms are within HB cut-off distance #CATENR_HB_CUTOFF.
 * @todo Since prefactor is given as macro it migh be removed from parameters ... unless same function being used for optimization purposes...
 */
double CATENR_HB_eq 		(double *r_H_i, double *r_N_i, double *r_O_j, double *r_C_j,double r,double pref);

/**
 * @brief Function calculate Caterpillar implicit solvatation energy
 *
 * Function substract number of contacts of given residue with other residues from treshold.
 * Treshold value being larger for hydrophobic residues, which if not burried (large number of contacts) would produce
 * energy penalty scaled by hydrophopbicity of residue *Hydrophob.
 */
double CATENR_Wat_eq		(int type, double contacts,double *Hydrophob,double *treshold,double pref);

/**
 * @brief Function calculate Caterpillar implicit solvatation energy of one residue
 *
 * Function substract number of contacts of given residue bond with other residues from treshold.
 * Treshold value being larger for hydrophobic residues, which if not burried (large number of contacts) would produce
 * energy penalty scaled by hydrophopbicity of residue *Hydrophob.
 * LT Tested OK 11/07/17 ( reproduces energy of a protein in input. Mind the prefactors!)
 *
 * @note Function iterate only on given protein chain, so that in case of multiple proteins
 * this will not calculate solvation energy of whole system but only within one protein!
 *
 * @todo Remove not used value n.
 */
double CATENR_Wat				(cat_prot *protein, int i, int n, double *Hydrophob,double *treshold,double pref);

/**
 * @brief Function calculate Caterpillar implicit solvatation energy of all residue
 *
 * @note Function iterate only on given protein chain, so that in case of multiple proteins
 * this will not calculate solvation energy of whole system but only within one protein!
 */
double CATENR_Wat_total		(cat_prot *protein, double *Hydrophob,double *treshold,double pref);

/**
 * @brief Function calculate bending energy
 */
double CATENR_Ca_Ca_Bend 	(double *r_Ca_i, double *r_Ca_j, double ref_bond, double pref);

/**
 * @brief Function read energy parameters from text file
 *
 * Function read Ca-Ca attractive interaction, hydrophobic scale, and hydrophobic treshold.
 *
 * LT Tested OK 08/05/17 (read wpot file correctly)
 */
void 	 CATENR_read_energy_params(int n, double **ptr_Hydrophob, double **ptr_H_treshold, double ***ptr_CaCa_Matrix, char filename[1024]);
#endif //__CAT_ENER
