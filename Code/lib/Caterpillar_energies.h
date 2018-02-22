#ifndef __CAT_ENER__
#define __CAT_ENER__

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include "Caterpillar.h"

#define CATENR_SAW_radius   4.0000000
#define CATENR_CaCa_RANGE   12.0000000
#define CATENR_CaCa_CUTOFF  17.0000000
#define CATENR_HB_CUTOFF		6.0000000
#define CATENR_HB_PREF   	  13.888	//3.1*8*0.56, with truncation er.
#define CATENR_CaCa_PREF    0.4400000000000000 	//1-0.56, with truncation er.
#define CATENR_Hydropathy   0.014655    					//Hydropathy scale from aapot
#define CATENR_WAT_PREF     0.00644819999999999926//CaCa_PREF*Hydropathy

#define CATENR_INFTY 100000000000
#define CATENR_ZERO_EN 0


double CATENR_bond 			(double r);
double CATENR_Saw_eq		(double r);
double CATENR_Ca_Ca_eq	(int type_i, int type_j, double bond, double **M, double pref);
double CATENR_Wat_eq		(int type, double contacts,double *Hydrophob,double *treshold,double pref);
double CATENR_Saw 			(double *r_Ca_i, double *r_Ca_j);
double CATENR_Saw_g 	    (double *r_atom_i, double *r_atom_j, double saw);
double CATENR_Saw_eq_g 	(double r, double saw);
double CATENR_Ca_Ca 		(double *r_Ca_i, int type_i, double *r_Ca_j, int type_j, int n, double **M,double pref);
double CATENR_HB  			(double *restrict r_H_i, double *restrict r_N_i, double *restrict r_O_j, double *restrict r_C_j, double pref);
double CATENR_HB_eq 		(double *r_H_i, double *r_N_i, double *r_O_j, double *r_C_j,double r,double pref);
double CATENR_Wat				(cat_prot *protein, int i, int n, double *Hydrophob,double *treshold,double pref);


double CATENR_Wat_total		(cat_prot *protein, double *Hydrophob,double *treshold,double pref);
double CATENR_Ca_Ca_Bend 	(double *r_Ca_i, double *r_Ca_j, double ref_bond, double pref);
void 	 CATENR_read_energy_params(int n, double **ptr_Hydrophob, double **ptr_H_treshold, double ***ptr_CaCa_Matrix, char filename[1024]);
#endif //__CAT_ENER
