#ifndef __ABM_PDBIO
#define __ABM_PDBIO
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//Return codes
#define PDBIO_SUCCESS 					(0)
#define PDBIO_ERROR_EMPTY_LINE	(1)
#define PDBIO_ERROR_SHORT_LINE	(2)
#define PDBIO_ERROR_NULL_PTR  	(3)
#define PDBIO_ERROR_OVERWRITE  	(4)
#define PDBIO_ERROR_WRITEMODE  	(5)

typedef struct {
	int	serial;
	char 			name[5];
	char 			altLoc;
	char 			resName[4];
	char 			chainID;
	int 	resSeq;
	char 			iCode;
	float 		x;
	float 		y;
	float 		z;
	float 		occupancy;
	float 		tempFactor;
	char 			element[3];
	char 			charge[3];
} pdb_atom;

struct _pdb_atom_list{
	pdb_atom values;
	struct _pdb_atom_list *prev_atom;
	struct _pdb_atom_list *next_atom;
	//struct _pdb_atom_list *last_atom;
};
typedef struct _pdb_atom_list pdb_atom_list;

/*put data on top of the stack*/
void PDBIO_atomlist_push 		( pdb_atom_list **head,  pdb_atom el );
/*recover data from top of the stack*/
void PDBIO_atomlist_pop 		( pdb_atom_list **head,  pdb_atom *el);

/* append a new element at the end of the list*/
void PDBIO_atomlist_append 	( pdb_atom_list **tail,  pdb_atom el );
/*delete entry from the list */
void PDBIO_atomlist_remove	( pdb_atom_list **prev,  pdb_atom_list *obj);


int PDBIO_read_atom  ( char *line , pdb_atom *at);
int PDBIO_write_atom ( pdb_atom *at, char *line );
int PDBIO_read_all_atoms (FILE * pdb_in, pdb_atom_list ** atoms);
int PDBIO_write_all_atoms ( FILE *pdb_out, pdb_atom_list *atoms);
//get things from structs
int PDBIO_atomlist_get_coord		 ( pdb_atom_list *atoms, int *N_at, double	 	 ***coord);
int PDBIO_atomlist_get_serials	 ( pdb_atom_list *atoms, int *N_at, int 			**serials);
int PDBIO_atomlist_get_names		 ( pdb_atom_list *atoms, int *N_at, char  		 ***names);
int PDBIO_atomlist_get_resNames	 ( pdb_atom_list *atoms, int *N_at, char 	***resNames);
int PDBIO_atomlist_get_altLoc		 ( pdb_atom_list *atoms, int *N_at, char 		**altLocs);
int PDBIO_atomlist_get_resSeq 	 ( pdb_atom_list *atoms, int *N_at, int  		**resSeqs);
int PDBIO_atomlist_get_occupancy ( pdb_atom_list *atoms, int *N_at, double		**occupancies);
int PDBIO_atomlist_get_tempFactor( pdb_atom_list *atoms, int *N_at, double		**tempFactors);
//
int PDBIO_coord2atomlist (double **coord, int *serial,char **name, char *altLoc, char **resName, int *resSeq, int N_at, pdb_atom chain_prop, pdb_atom_list **atoms);
//Utils
void strip_spaces    ( char *s);

#endif // __ABM_PDBIO
