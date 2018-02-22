#ifndef __CAT_IO
#define __CAT_IO
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ABM_pdbio.h"
#include "Caterpillar.h"


void CATIO_fastadecoder ( int *Dec, char *Enc, int Seq_Length);
void CATIO_cat2pdb 			( char *filename, char *writemode, char *remark, cat_prot *prot, int N_prot );
//void CATIO_pdb2cat_skip_CB  		( cat_prot **prot, int *N_prot, char *filename);
void CATIO_pdb2cat_skip_CB  		( cat_prot **prot, int *N_prot, char *filename);
void CATIO_pdb2cat_keep_CB  		( cat_prot **prot, int *N_prot, char *filename);
//void CATIO_pdb2cat  	( cat_prot **prot, int *N_prot, int N_atom_types, char ATnames[N_atom_types][5], char Loc_keep, char *filename);
void CATIO_pdb2cat  	( cat_prot **prot, int *N_prot, int N_atom_types, char ATnames[N_atom_types][5], char Loc_keep, char *filename);
void CATIO_cat2oldbin 		( char *filename, char *writemode, cat_prot *prot, int N_prot );
void CATIO_oldbin2cat  		( cat_prot *prot, int N_prot, char *filename);

#endif // __CAT_IO
