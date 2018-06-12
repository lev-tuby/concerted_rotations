#ifndef __ABM_PDBIO
#define __ABM_PDBIO
#include <stdio.h>
/**
 * @file
 * @brief Header file for all functions to work with PDB formated atom data
 * Contain function for reading and writing PDB files
 * @todo Functions for getting different parts of pdb_atom_list are redundant, substitute with a single function.
 * @todo pdb_atom list might contain number of atoms there ... 
 */

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
/** @brief Value returned if PDB line is shorther then format defined length */
#define PDBIO_ERROR_SHORT_LINE	(2)
/** @brief Value returned if some of the input parameters are NULL pointers */
#define PDBIO_ERROR_NULL_PTR  	(3)
/** @brief Value returned if pointer to output values contains data (would cause mem leak!) */
#define PDBIO_ERROR_OVERWRITE  	(4)
/** @brief Value not used */
#define PDBIO_ERROR_WRITEMODE  	(5)
/** @} */ // end of PDBIO_RETURN_CODES

/**
 * @brief Structure storing the PDB parameters of a single atom
 */
typedef struct {
    int     serial;         /**< Atom index. */
    char    name[5];        /**< Atom name. */
    char    altLoc;         /**< Alternative loacation. */
    char    resName[4];     /**< Residue name. */
    char    chainID;        /**< Chain ID. */
    int     resSeq;         /**< Residue index. */
    char    iCode;          /**< iCode. */
    // coordinates
    float   x;              /**< Atom x coordinate. */
    float   y;              /**< Atom y coordinate. */
    float   z;              /**< Atom z coordinate. */
    float   occupancy;      /**< Occupancy factor. */
    float   tempFactor;     /**< Temperature factor. */
    char    element[3];     /**< Element type. */
    char    charge[3];      /**< Atomic charge. */
} pdb_atom;

/**
 * @brief Structure storing multiple atoms in an ordered manner
 *
 * Doubly-linked list of pdb_atom structures.
 */
struct _pdb_atom_list{
	pdb_atom values;                    /**< Actual atom structure in the linked list. */
	struct _pdb_atom_list *prev_atom;   /**< Link to previous atom. */
	struct _pdb_atom_list *next_atom;   /**< Link to next atom. */
	//struct _pdb_atom_list *last_atom;
};
typedef struct _pdb_atom_list pdb_atom_list;

/**
 * @brief Test function: read data from a PDB and print it to a second one.
 *
 */
int PDBIO_test (const char *input_path, const char *output_path);

/**
 * @brief Remove spaces from beginning and end of a string
 */
void strip_spaces    (char *s);

/**
 * @brief Reads one line from a  pdb and convert it to a #pdb_atom structure
 */
int PDBIO_read_atom  (char *line, pdb_atom *at);

/**
 * @brief Writes one atom line
 *
 * Alignment of one-letter atom name such as C starts at column 14, while two-letter atom name such as FE starts at column 13.
 */
int PDBIO_write_atom (pdb_atom *at, char *line );

/**
 * @brief Reads all atoms from a pdb file and store them in a #_pdb_atom_list structure
 *
 */
int PDBIO_read_all_atoms (FILE * pdb_in, pdb_atom_list ** atoms);

/**
 * @brief Writes a pdb file starting from a #_pdb_atom_list structure
 *
 */
int PDBIO_write_all_atoms ( FILE *pdb_out, pdb_atom_list *atoms);

/**
 * @brief Pushes a new element on top of the list in the #_pdb_atom_list structure
 *
 */
void PDBIO_atomlist_push 		( pdb_atom_list **head,  pdb_atom el );

/**
 * @brief Pops the top element from the list in the #_pdb_atom_list structure
 *
 * @note UNUSED. NOT TESTED
 */
void PDBIO_atomlist_pop 		( pdb_atom_list **head,  pdb_atom *el);

/**
 * @brief Appends a new element to the end of the list in the #_pdb_atom_list structure
 *
 */
void PDBIO_atomlist_append 	( pdb_atom_list **tail,  pdb_atom el );

/**
 * @brief Removes a given element from the list in the #_pdb_atom_list structure
 * @note UNUSED. NOT TESTED
 *
 */
void PDBIO_atomlist_remove	( pdb_atom_list **prev,  pdb_atom_list *obj);


//get things from structs
/**
 * @brief Gets number and coordinates from a #_pdb_atom_list structure
 *
 */
int PDBIO_atomlist_get_coord		 ( pdb_atom_list *atoms, int *N_at, double	 	 ***coord);

/**
 * @brief Extracts an array of atom indexes from #_pdb_atom_list
 *
 * Iterates through a #_pdb_atom_list and read array of int values of atom indexes for each atom.
 * Does not reorder values in any way.
 */
int PDBIO_atomlist_get_serials	 ( pdb_atom_list *atoms, int *N_at, int 			**serials);

/**
 * @brief Extracts an array of atom names from a #_pdb_atom_list
 *
 * Iterates through a #_pdb_atom_list and read array of *char of atom names for each atom.
 * Does not reorder values in any way.
 */
int PDBIO_atomlist_get_names		 ( pdb_atom_list *atoms, int *N_at, char  		 ***names);

/**
 * @brief Extracts an array of residue names from a #_pdb_atom_list
 *
 * Iterates through a #_pdb_atom_list and read array of *char of residue names for <b>each atom</b>.
 * Does not reorder values in any way.
 */
int PDBIO_atomlist_get_resNames	 ( pdb_atom_list *atoms, int *N_at, char 	***resNames);

/**
 * @brief Extracts an array of altLocs from a #_pdb_atom_list
 *
 * Iterates through a #_pdb_atom_list and read array of char of altLocs for <b>each atom</b>.
 * Does not reorder values in any way.
 */
int PDBIO_atomlist_get_altLoc		 ( pdb_atom_list *atoms, int *N_at, char 		**altLocs);

/**
 * @brief Extracts an array of residue numbers from a #_pdb_atom_list
 *
 * Iterates through a #_pdb_atom_list and read array of int of residue numbers for <b>each atom</b>.
 * Does not reorder values in any way.
 * Tested OK. LT 01.08.16
 */
int PDBIO_atomlist_get_resSeq 	 ( pdb_atom_list *atoms, int *N_at, int  		**resSeqs);

/**
 * @brief Extracts an array of occupancies from a #_pdb_atom_list
 *
 * Iterates through a #_pdb_atom_list and read array of doubles of occupancies for <b>each atom</b>.
 * Does not reorder values in any way.
 */
int PDBIO_atomlist_get_occupancy ( pdb_atom_list *atoms, int *N_at, double		**occupancies);

/**
 * @brief Extracts an array of temperature factors from a #_pdb_atom_list
 *
 * Iterates through a #_pdb_atom_list and read array of doubles of temperature factors for <b>each atom</b>.
 * Does not reorder values in any way.
 */
int PDBIO_atomlist_get_tempFactor( pdb_atom_list *atoms, int *N_at, double		**tempFactors);

/**
 * @brief Takes an ordered array of atom parameters and fill them from a #_pdb_atom_list structure
 *
 * Takes separate arrays for each column in pdb (coordinates, atom index, atom name, residue name etc.) and fill them in a #_pdb_atom_list.
 * If data for alternative location are not supplied (*altLoc==NULL), then empty space is leaved in pdb.
 * If data for atom indices are not supplied (*serial==NULL), then atoms are numbered from 1.
 *
 * @note parameters in all <b>arrays have to be ordered</b>.
 */
int PDBIO_coord2atomlist (double **coord, int *serial,char **name, char *altLoc, char **resName, int *resSeq, int N_at, pdb_atom chain_prop, pdb_atom_list **atoms);


#endif // __ABM_PDBIO
