#ifndef __CATM_H__
#define __CATM_H__
/**
 * @file
 * @brief Header file for moves on Catepillar model
 */

#include "Caterpillar.h"
#include <gsl/gsl_rng.h>
#include "../generic_move/CR_precomp.h"

/**
 * @brief Stores the indeces of residue pairs whose distances were modified by a move
 *
 * Used to speed up energy calculation.
 * @note #dist is actually always set to -1.0 and not used in the current implementation
 */
typedef struct {
    int i_1;        /**< Index of first residue. */
    int i_2;        /**< Index of second residue. */
    double dist;    /**< Unused */
} resid_pair;

/**
 * @brief Structure to store data on MC moves
 */
typedef struct {
	int _l1;                /**< Maximum number of residues that can be modified by a single move. */
	int _l2;                /**< Maximum number of residues on the backbone. */
	int N_moved;            /**< Number of residues modified by the last move. */
	int *moved_res;         /**< Array storing the indeces of residues that were moved in the last move. */
	int N_pairs;            /**< Number of residues pairs whose distances where modified by the last move (#mod_pairs ' dimension ).*/
	resid_pair *mod_pairs;   /**< Array storing pairs of residues whose distances were modified by the last move. */
} mc_move_data;

/**
 * @brief Allocates a #mc_move_data structure
 */
mc_move_data * CATMV_mc_move_data_alloc(int max_move_size,int len);

/**
 * @brief Deallocates a #mc_move_data structure
 */
void CATMV_mc_move_data_free( mc_move_data *mvdt);

/**
 * @brief Performs a crankshaft rotation of part of a #cat_prot backbone
 */
void CATMV_cranck	(mc_move_data *cranck_data, cat_prot * p, gsl_rng *rng_r);

/**
 * @brief Performs a pivot rotation of part of a #cat_prot backbone
 */
void CATMV_pivot	(mc_move_data *pivot_data, 	cat_prot *p, 	gsl_rng *rng_r, int max_move_size);

/**
 * @brief Performs a concerted rotation of 4 consecutive a.a. of a #cat_prot backbone
 */
void CATMV_concerted_rot(mc_move_data *ra_data, cat_prot *p, gsl_rng * rng_r, double sigma);

/**
 * @brief performs a random concerted rotation of a set of atoms
 * @note to be moved back to the .c file asap. Placed here only for testing purposes.
 */
int random_rot(cr_input_data bb_out, cr_input_data bb_in, gsl_rng *rng_r, double sigma);

#endif //__CATM_H__
