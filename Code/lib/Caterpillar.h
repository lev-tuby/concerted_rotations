#ifndef __CATERPILLAR__
#define __CATERPILLAR__
/**
 * @file
 * @brief Function for manipulation of Cterpillar peptide backbone
 * @todo Remove deprecated functions and decide which functions will stay or which functions will split etc.
 */

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>


/** @brief Size of alphabet */
#define CAT_S 21

/** @brief Acceptable residues names
 *
 * DUM for Dummy residue.
 *
*/
#define CAT_AADefs {"ALA","CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE",\
	"LYS" ,"LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL",\
	"TRP", "TYR","DUM"}

/** @brief Ordered string of one char FASTA code for aminoacids
 *
 * Residues in the code are stored as int values. Value that coresponds to given one letter code is given in #CAT_AACODE, where number coresponding to n-th one letter code is CAT_AACODE[n].
*/
#define CAT_FASTA "ACDEFGHIKLMNPQRSTVWYX" 

/** @brief Array that define decoding of amino acids to number indexes
 *
 * Residue types are represented in code by simple int values, where to each one leter FASTA code one number coresponds. Order of fasta codes can be deduced from #CAT_FASTA.
 * 
 * This differ from ivan's fasta decoder. The reason is that it makes it much simpler to assign the code to the letter and viceversa (the num. id. is simplythe index of the letter in the #CAT_FASTA string).
 *
*/
#define CAT_AACODE {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20}

/** @defgroup PDBIO_RETURN_CODES
 *  Indices (order) of atoms in cat_prot structure
 *  @{
 */
/** @brief Index of hydrogen */
#define ATOM_H 0
/** @brief Index of Nitrogen */
#define ATOM_N 1
/** @brief Index of Ca carbon */
#define ATOM_CA 2
/** @brief Index of C */
#define ATOM_C 3
/** @brief Index of oxygen */
#define ATOM_O 4
/** @brief Index of Cb carbon */
#define ATOM_CB 5
/** @} */ // end of PDBIO_RETURN_CODES

//Atom names
/** @brief Atom names used in cod
 *  @Note Order is not same as one used in code ... in code first atom is N
 */
#define CAT_ATnames {"H","N","CA","C","O","CB"}

/** @brief Number of atoms per residue */
#define NATOM 6

/** @brief By which we move along nromal to C-CA-N plane in construction of CB atoms */
#define CAT_Cb_AZIMUTH 1.2

/** @brief Length of Ca-Cb bond
 * @todo Change its name to CAT_Rbond_CaCb ?
 */
#define CAT_Rbond_CCb 1.000000

/** @defgroup Peptide backbone geometric parameters
 *  Bond lenghts, bond angles etc.
 *  @{
 */

/** @brief H-C-N andgle dot product */
#define CAT_Plat__C_HdotC_N (2.4318886938725948)

/** @brief C-N-Ca andgle dot product */
#define CAT_Plat__N_CAdotC_N  (2.8923429607978419)

/** @brief O-C-H andgle dot product */
#define CAT_Plat__C_HdotC_O  (-2.1259814163093753)

/** @brief Ca-N bond length
 * sist->Rbond[0]= 1.45000
*/
#define CAT_Rbond_CaN 1.4500000

/** @brief C-Ca bond length
 * sist->Rbond[1]=1.52000
*/
#define CAT_Rbond_CCa 1.5200000

/** @brief C-O bond length
 * sist->Rbond[2]= 1.23000
*/
#define CAT_Rbond_CO  1.2300000

/** @brief C-N bond length
 * sist->Rbond[3]= 1.33000
*/
#define CAT_Rbond_CN  1.3300000

/** @brief N-H bond length
 * sist->Rbond[4]= 1.00000
*/
#define CAT_Rbond_NH  1.0000000

/** @brief N-C bond length
 * sist->Rbond[5], used for Bending instead of the angle..
*/
#define CAT_Rbond_NC  2.4479800


/** @brief N-Ca-C angle */
#define CAT_angle_NCaC 1.9373154697137058

/** @brief Ca-C-N angle */
#define CAT_angle_CaCN 2.017600615305445

/** @brief C-N-Ca angle */
#define CAT_angle_CNCa 2.1275563581810877

/** @brief C-N-H angle */
#define CAT_angle_CNH  2.0856684561332237

/** @brief Ca-C-O angle */
#define CAT_angle_CaCO 2.1135937241651326
/** @} */ // end of Peptide backbone geometric parameters

/**
 * @brief Main structure to store protein configuration
 *
 * Data structure store coordinates of all atoms in protein, number of residues, torsion angles.
 * Data stracture is used for energy calculation and is altered in moves.
 */
typedef struct {
	size_t n_atom_per_res;  /**< Number of atoms per residue (6 if CB included otherwise 5). */
	size_t n_atoms;         /**< Number of all atoms in protein. */
	size_t n_res;           /**< Number of all residues in protein. */
	double **coord;         /**< Array of coordinates of all atoms in protein. */
	double ** N;            /**< Array slice of coord array where all positions of N atoms are stored. N[i] corespond to N atom in i-th residue. */
	double **CA;            /**< Array slice of coord array where all positions of Ca atoms are stored. Ca[i] corespond to Ca atom in i-th residue. */
	double ** C;            /**< Array slice of coord array where all positions of C atoms are stored. C[i] corespond to C atom in i-th residue. */
	double ** O;            /**< Array slice of coord array where all positions of O atoms are stored. O[i] corespond to O atom in i-th residue. */
	double ** H;            /**< Array slice of coord array where all positions of H atoms are stored. H[i] corespond to H atom in i-th residue. */
	double **CB;            /**< Array slice of coord array where all positions of Cb atoms are stored. Cb[i] corespond to Cb atom in i-th residue. */
	double *phi;            /**< Array of all phi thorsion angles. phi[i] corespond to phi dihedral in i-th residue. */
	double *psi;            /**< Array of all psi thorsion angles. psi[i] corespond to psi dihedral in i-th residue. */
	double *contacts;       /**< Number of contacts for each residue. Value is used for wanter contacts energy. */
	int    *residues;       /**< Array contain residue types of all residues. Values are determined by #CAT_AACODE where translation to nonce letter FASTA code is given in #CAT_FASTA. */
} cat_prot;

/**
 * @brief Function initialize cat_prot strucure
 *
 * All coordinates are set to 0.0, all dihedrals are set to 0.0. p->=n_res*n_atom_per_res, p->n_res=n_res, p->n_atom_per_res=n_atom_per_res, p->residues[i]=-1, p->contacts[i]=-1.0.
 * Atoms are aliased ... it is kind of dangerous though. It remains unclear how optimization is affected inside a function.
 */
cat_prot * CAT_prot_alloc (size_t n_res, size_t n_atom_per_res);

/**
 * @brief Deallocation of cat_prot structure
 */
void CAT_prot_free (cat_prot * pr);

/**
 * @brief Function print in stdout position of all atoms in protein
 */
void CAT_print(cat_prot *protein);

/**
 * @brief Function build up cat_prot from dihedral angles
 *
 * Build a Caterpillar protein, assuming that the first dihedral angle passed through dihed is a phi angle.
 */
cat_prot * CAT_build_from_dihed ( int n_res, size_t n_atom_per_res, double *orig, double *dihed, char *seq);

/**
 * @brief Function build up liner protein chain
 *
 * Tested OK. LT 02.08.16
 */
void CAT_set_prot_linear ( cat_prot * protein, double *orig, double alpha );
//void CAT_set_residues ( cat_prot * protein, int *Dec);

/**
 * @brief Function set protein residue types from FASTA string
 */
void CAT_set_residues_fasta ( cat_prot * protein, int Seq_Length, char *Enc);

/**
 * @brief Function add hydrogens to existing peptide backbone
 */
void CAT_insert_hydrogens ( cat_prot * protein );

/**
 * @brief Function add pozition o CB atom to given residue of protein
 */
void CAT_insert_cbeta ( cat_prot *protein, int res_i, double azimut, double bond_length	);

/**
 * @brief Function correct backbone distortions
 *
 * Tested OK. 03.08.16. (A shrinked Caterpillar protein is reinflated to itself)
 */
void CAT_rescale ( cat_prot * protein );
//void CAT_add_peptide ( cat_prot * protein, int i, double phi, double psi);
//void CAT_add_peptide ( cat_prot * protein, int I, double phi, double psi, double angle_NCaC);

/**
 * @brief Function modifie residue configuration base on dihedrals
 *
 * Function reconstruct given residue based on dihedral angles in cat_prot structure.
 * Tested OK. LT 02.08.16
 */
int CAT_add_peptide ( cat_prot *p, int I, double phi, double alpha, double psi );

/**
 * @brief Function recalculate protein dihedral angles
 *
 * Function iterate through whole protein and recalculate p->phi and p->psi.
 */
void CAT_prot_dihedrals (cat_prot *protein);

/**
 * @brief Function calculate phi angle of residue in protein
 */
double compute_phi(cat_prot *p,int c);

/**
 * @brief Function calculate psi angle of residue in protein
 */
double compute_psi(cat_prot *p,int c);

/**
 * @brief Function recalculate protein dihedral angles
 *
 * Function iterate through whole protein and recalculate p->phi and p->psi.
 * @todo have to be rewritten
 *
 * @param[in,out]   *p        Protein which dihedrals are recalculated
 *
 * @return \c void
 */
void compute_dihedrals(cat_prot *p);

/**
 * @brief Function calculate dihedral from 4 atoms
 *
 * This function was taken from LAMMPS.
 */
double calc_dihedralf_angle(double *atom_1, double *atom_2, double *atom_3, double *atom_4);

/**
 * @brief Don not know ... Luca?
 */
void build_peptide ( gsl_matrix *pep);
#endif //__CATERPILLAR__
