#ifndef __CATM_H__
#define __CATM_H__
/**
 * @file
 * @brief Header file for moves on Catepillar model
 */

#include <stdio.h>
#include <math.h>
#include "Caterpillar.h"
#include <gsl/gsl_rng.h>
#include "../generic_move/CR_precomp.h"

/**
 * @brief Pairs of atoms which distances were changed
 *
 * Pair contain indexes of atoms which distances have changed and their current distance. Help in interaction energy calculation.
 * @todo Mayber rename to residue pairs ... since they are used as residue pairs and #dist is actualy not used anyway ...
 * @note #dist is actually always set to -1.0 and not used
 */
typedef struct {
    int i_1;        /**< Index of first atom. */
    int i_2;        /**< Index of second atom. */
    double dist;    /**< Current distance between atoms. */
} atom_pair;

/**
 * @brief Structure to store data on MC moves
 */
typedef struct {
	int _l1;                /**< Maximal number of residues to bee moved during one move. */
	int _l2;                /**< Maximal number of residues in sequnce. */
	int N_moved;            /**< Number of residues moved in in the last move. */
	int *moved_res;         /**< Array of indexes of residues that were moved in the last move. */
	int N_pairs;            /**< Number of residues that were moved in the last move (only data that are lower than this number are relevant to last move!). */
	atom_pair *mod_pairs;   /**< Array of pairs of rasidues which distances changed in the last move. */
} mc_move_data;

/**
 * @brief Function allocate #mc_move_data structure
 *
 * #atom_pair that are moved are initialized to {-1, -1, -1.0}.
 *
 * @note moved_res are not initialized
 */
mc_move_data * CATMV_mc_move_data_alloc(int max_move_size,int len);

/**
 * @brief Function deallocate #mc_move_data structure
 */
void CATMV_mc_move_data_free( mc_move_data *mvdt);

/**
 * @brief Function perform crankshaft move on protein
 */
void CATMV_cranck	(mc_move_data *cranck_data, cat_prot * p, gsl_rng *rng_r);

/**
 * @brief Function perform pivot move on protein
 *
 * Pivot move here operate only on ends of protein. How many residues from both ends are affected is determined by #max_move_size.
 *
 * @note Tested OK with VMD
 */
void CATMV_pivot	(mc_move_data *pivot_data, 	cat_prot *p, 	gsl_rng *rng_r, int max_move_size);

/**
 * @brief Function perform concerted rotation move on protein
 */
void CATMV_concerted_rot(mc_move_data *ra_data, cat_prot *p, gsl_rng * rng_r, double sigma);

//Move those back to the .c file asap is there due to need of testing Concerted rotation ...
int random_rot(cr_input_data bb_out, cr_input_data bb_in, gsl_rng *rng_r, double sigma);

#endif //__CATM_H__
