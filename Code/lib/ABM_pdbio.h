#ifndef __ABM_PDBIO
#define __ABM_PDBIO
/**
 * @file
 * @brief Header file for all functions to work with PDB formated atom data
 * Contain function for reading and writing PDB files
 * @todo Functions for getting diferent parts of pdb_atom_list ... are well redundant there might be one function that does that all ...
 * @todo Also ... pdb_atom list might contain number of atoms there ... 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/** @defgroup PDBIO_RETURN_CODES
 *  All posible return codes from PDBIO functions
 *  @{
 */
/** @brief Default value used in intialization of PDBIO return values */
#define PDBIO_DEFAULT           (-1)
/** @brief Value returned if function succeeded */
#define PDBIO_SUCCESS 			(0)
/** @brief Value returned if PDB line is empty */
#define PDBIO_ERROR_EMPTY_LINE	(1)
/** @brief Value returned if PDB line is shorther then formate defined length */
#define PDBIO_ERROR_SHORT_LINE	(2)
/** @brief Value returned if some of input parameters are NULL pointer */
#define PDBIO_ERROR_NULL_PTR  	(3)
/** @brief Value returned if pointer to output values contain data (would cause mem leak!) */
#define PDBIO_ERROR_OVERWRITE  	(4)
/** @brief Value not used ... ??? */
#define PDBIO_ERROR_WRITEMODE  	(5)
/** @} */ // end of PDBIO_RETURN_CODES

/**
 * @brief Structure for holding PDB parameters of single atom
 */
typedef struct {
    int     serial;
    char    name[5];
    char    altLoc;
    char    resName[4];
    char    chainID;
    int     resSeq;
    char    iCode;
    // coordinates
    float   x;
    float   y;
    float   z;
    float   occupancy;
    float   tempFactor;
    char    element[3];
    char    charge[3];
} pdb_atom;

/**
 * @brief Structure for holding multiple atoms in order manner
 *
 * Simple linked list of pdb_atom structures.
 */
struct _pdb_atom_list{
	pdb_atom values;
	struct _pdb_atom_list *prev_atom;
	struct _pdb_atom_list *next_atom;
	//struct _pdb_atom_list *last_atom;
};
typedef struct _pdb_atom_list pdb_atom_list;

/**
 * @brief Function test functions in the PDBIO lib
 *
 * Function test functions by reading data from PDB file and printing them out to other pdb file.
 */
int PDBIO_test (const char *input_path, const char *output_path);

/**
 * @brief Function for removing space at the end of string
 *
 * Quite clever trimming loop I've found at http://stackoverflow.com/questions/7775138/strip-whitespace-from-a-string-in-place
 * The termination condition tests that s[i]==NULL and set s[j]=s[i] in one swoop.
 * Tested OK. LT 01.08.16
 */
void strip_spaces    (char *s);

/**
 * @brief Function read one line of pdb and convert it to pdb_atom structure
 *
 * ...
 * Tested OK. LT 02.08.16
 */
int PDBIO_read_atom  (char *line, pdb_atom *at);

/**
 * @brief Function read one line of pdb and convert it to pdb_atom structure
 *
 * From http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOMge
 * Alignment of one-letter atom name such as C starts at column 14, while two-letter atom name such as FE starts at column 13.
 * Tested OK. LT 01.08.16
 */
int PDBIO_write_atom (pdb_atom *at, char *line );

/**
 * @brief Function take FILE structure of pdb and read it in pdb_atom_list structure
 *
 * Tested OK. LT 01.08.16 -- NOT VALID ANYMORE
 */
int PDBIO_read_all_atoms (FILE * pdb_in, pdb_atom_list ** atoms);

/**
 * @brief Function take FILE structure of output pdb and print content of pdb_atom_list structure
 *
 * Tested OK. LT 01.08.16
 */
int PDBIO_write_all_atoms ( FILE *pdb_out, pdb_atom_list *atoms);

/**
 * @brief Function push a new element on top of the list in the pdb_atom_list structure
 *
 * NOT TESTED
 */
void PDBIO_atomlist_push 		( pdb_atom_list **head,  pdb_atom el );

/**
 * @brief Function recover top element from the list in the pdb_atom_list structure
 *
 * NOT TESTED
 */
void PDBIO_atomlist_pop 		( pdb_atom_list **head,  pdb_atom *el);

/**
 * @brief Function append new element to the end of the list in pdb_atom_list structure
 *
 * NOT TESTED
 */
void PDBIO_atomlist_append 	( pdb_atom_list **tail,  pdb_atom el );

/**
 * @brief Function remove element from the list in pdb_atom_list structure
 *
 * NOT TESTED
 */
void PDBIO_atomlist_remove	( pdb_atom_list **prev,  pdb_atom_list *obj);


//get things from structs
/**
 * @brief Function return number and coordinates from the pdb_atom_list structure
 *
 * Function iterate through whole pdb_atom_list and count and store PDB coordinates in d2t(N, 3) array
 * Note: should be able to doo whole stuff in one for cycle ... also i variable is not used ...
 * Tested OK. LT 01.08.16
 */
int PDBIO_atomlist_get_coord		 ( pdb_atom_list *atoms, int *N_at, double	 	 ***coord);

/**
 * @brief Function return ordered array of atom indexes in pdb_atom_list
 *
 * Function iterate through whole pdb_atom_list and read array of int values of atom indexes for each atom in pdb_atom_list.
 * Function do not reorder values in any way.
 * Tested OK. LT 01.08.16
 */
int PDBIO_atomlist_get_serials	 ( pdb_atom_list *atoms, int *N_at, int 			**serials);

/**
 * @brief Function return ordered array of atom names in pdb_atom_list
 *
 * Function iterate through whole pdb_atom_list and read array of *char of atom names for each atom in pdb_atom_list.
 * Function do not reorder values in any way.
 * Tested OK. LT 01.08.16
 */
int PDBIO_atomlist_get_names		 ( pdb_atom_list *atoms, int *N_at, char  		 ***names);

/**
 * @brief Function return ordered array of residue names in pdb_atom_list
 *
 * Function iterate through whole pdb_atom_list and read array of *char of residue names for <b>each atom</b> in pdb_atom_list.
 * Function do not reorder values in any way.
 * Tested OK. LT 01.08.16
 */
int PDBIO_atomlist_get_resNames	 ( pdb_atom_list *atoms, int *N_at, char 	***resNames);

/**
 * @brief Function return ordered array of altLocs in pdb_atom_list
 *
 * Function iterate through whole pdb_atom_list and read array of char of altLocs for <b>each atom</b> in pdb_atom_list.
 * Function do not reorder values in any way.
 * Tested OK. LT 01.08.16
 */
int PDBIO_atomlist_get_altLoc		 ( pdb_atom_list *atoms, int *N_at, char 		**altLocs);

/**
 * @brief Function return ordered array of residue numbers in pdb_atom_list
 *
 * Function iterate through whole pdb_atom_list and read array of int of residue numbers for <b>each atom</b> in pdb_atom_list.
 * Function do not reorder values in any way.
 * Tested OK. LT 01.08.16
 */
int PDBIO_atomlist_get_resSeq 	 ( pdb_atom_list *atoms, int *N_at, int  		**resSeqs);

/**
 * @brief Function return ordered array of ocupancies in pdb_atom_list
 *
 * Function iterate through whole pdb_atom_list and read array of doubles of ocupancies for <b>each atom</b> in pdb_atom_list.
 * Function do not reorder values in any way.
 */
int PDBIO_atomlist_get_occupancy ( pdb_atom_list *atoms, int *N_at, double		**occupancies);

/**
 * @brief Function return ordered array of temperature factors in pdb_atom_list
 *
 * Function iterate through whole pdb_atom_list and read array of doubles of temperature factors for <b>each atom</b> in pdb_atom_list.
 * Function do not reorder values in any way.
 */
int PDBIO_atomlist_get_tempFactor( pdb_atom_list *atoms, int *N_at, double		**tempFactors);

/**
 * @brief Function take ordered array of atom parameters and fill them in pdb_atom_list structure
 *
 * Function take separate arraies for each column in pdb (coordinates, atom index, atom name, residue name etc.) and fill them in pdb_atom_list.
 * If data for alternative location are not supplied (*altLoc==NULL), then empty space is leaved in pdb.
 * If data for atom indices are not supplied (*serial==NULL), then atoms are numbered from 1.
 * Tested OK. LT 01.08.16
 *
 * Note: parameters in all <b>arraies have to be ordered</b>.
 * TODO: Only existence of coord is tested but not of other input arrays.
 *       What is the purpose of chain_prop? it is nor really used ...
 */
int PDBIO_coord2atomlist (double **coord, int *serial,char **name, char *altLoc, char **resName, int *resSeq, int N_at, pdb_atom chain_prop, pdb_atom_list **atoms);


#endif // __ABM_PDBIO
