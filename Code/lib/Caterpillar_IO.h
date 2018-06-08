#ifndef __CAT_IO__
#define __CAT_IO__
/**
 * @file
 * @brief Header file for reading/writting PDB files
 *
 * @note Most of functions are there for compatibility reasons.
 *
 * @todo what about having one function that save/read configuration to/from PDB and only macros representing format are used ... would hide internal layer
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ABM_pdbio.h"
#include "Caterpillar.h"

/**
 * @brief Function code FASTA code in number code used internaly in Catepillar model
 *
 * @note Potentionaly dangerous function which might cause writing in non-allocated memory etc.
 *
 * @todo I would remove function from here since now coding/decoding is done with #CAT_FASTA and #CAT_AACODE.
 */
void CATIO_fastadecoder ( int *Dec, char *Enc, int Seq_Length);

/**
 * @brief Function take internal reresentation of Caterpillar backbone and save it in PDB file
 *
 * @todo ... change usage of N_prot ... its not used corectly.
 */
void CATIO_cat2pdb 			( char *filename, char *writemode, char *remark, cat_prot *prot, int N_prot );

/**
 * @brief Function read PDB file and save it to internal Caterpillar model
 *
 * Tested OK 12/07/17
 * @note Does not insert hydrogens, does not rescale, does not compute dihedrals
 * @todo ... change usage of N_prot ... is not usd in function at all.
 */
void CATIO_pdb2cat  	( cat_prot **prot, int *N_prot, int N_atom_types, char ATnames[N_atom_types][5], char Loc_keep, char *filename);

/**
 * @brief Function read PDB file and save it to internal Caterpillar model (keep CB atoms positions)
 *
 * @todo ... change usage of N_prot ... is not usd in function at all.
 */
void CATIO_pdb2cat_keep_CB  		( cat_prot **prot, int *N_prot, char *filename);

/**
 * @brief Function take internal reresentation of Caterpillar backbone and save it in binary format of PDB file
 *
 * @todo ... change usage of N_prot ... its not used corectly.
 */
void CATIO_cat2oldbin 		( char *filename, char *writemode, cat_prot *prot, int N_prot );

/**
 * @brief Function read binary PDB file format and save it to internal Caterpillar model
 *
 * Compatibility version. N_prot is passed to the function. It would make
 * much more sense to save it on the binary file. Also prot must already be
 * allocate with a valid number of atoms and residues, which will be read from
 * the file.
 */
void CATIO_oldbin2cat  		( cat_prot *prot, int N_prot, char *filename);

#endif // __CAT_IO__
