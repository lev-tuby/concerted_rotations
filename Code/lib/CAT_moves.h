#ifndef __CATM_H
#define __CATM_H

#include <stdio.h>
#include <math.h>
#include "Caterpillar.h"
#include <gsl/gsl_rng.h>
#include "../generic_move/CR_precomp.h"


typedef struct {
	int i_1;
	int i_2;
	double dist;
} atom_pair;

typedef struct {
	int _l1;
	int _l2;
	int N_moved;
	int *moved_res;
	int N_pairs;
	atom_pair *mod_pairs;
}mc_move_data;

mc_move_data * CATMV_mc_move_data_alloc(int max_move_size,int len);
void CATMV_mc_move_data_free( mc_move_data *mvdt);
void CATMV_cranck	(mc_move_data *cranck_data, cat_prot * p, gsl_rng *rng_r);
void CATMV_pivot	(mc_move_data *pivot_data, 	cat_prot *p, 	gsl_rng *rng_r, int max_move_size);
void CATMV_concerted_rot(mc_move_data *ra_data, cat_prot *p, gsl_rng * rng_r, double sigma);
//Move those back to the .c file asap
int random_rot(cr_input_data bb_out, cr_input_data bb_in, gsl_rng *rng_r, double sigma);

#endif //__CATM_H
