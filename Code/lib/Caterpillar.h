#ifndef __CATERPILLAR__
#define __CATERPILLAR__
/**
 * @file
 * @brief Interface for the Caterpillar protein model. 
 * @todo Remove deprecated functions and decide which functions will stay or which functions will split etc.
 */

#include <gsl/gsl_matrix.h>


/** @brief alphabet Size (number of "amino-acid" types in the interaction matrix)*/
#define CAT_S 21

/** @brief Acceptable residues names
 *
 * DUM for Dummy residue.
 *
*/
#define CAT_AADefs {"ALA","CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE",\
	"LYS" ,"LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL",\
	"TRP", "TYR","DUM"}

/** @brief Ordered string of 1-letter FASTA codes for aminoacids types
 *
 * The corresponding int values used by the code to recognize different
 * aminoacids are given by the number found in CAT_AACODE at the corresponding 
 * index.
 *
*/
#define CAT_FASTA "ACDEFGHIKLMNPQRSTVWYX" 

/** @brief Array defining the internal integer values corresponding to each 
 * 1-letter code in CAT_FASTA.
 *
 * The simplest version is to have the values stored in CAT_AACODE correspond 
 * with the index: CAT_AACODE[j]=j.
 *
*/
#define CAT_AACODE {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20}

/** @defgroup PDBIO_RETURN_CODES
 *  Indices (order) of atoms in the caterpillar data structure `cat_prot`
 *  @{
 */
/** @brief Index for hydrogen */
#define ATOM_H 0
/** @brief Index for Nitrogen */
#define ATOM_N 1
/** @brief Index for Ca carbon */
#define ATOM_CA 2
/** @brief Index for C */
#define ATOM_C 3
/** @brief Index for oxygen */
#define ATOM_O 4
/** @brief Index for Cb carbon */
#define ATOM_CB 5
/** @} */ // end of PDBIO_RETURN_CODES

//Atom names
/** @brief Atom names used in code
 *  @Note Order is not same as one used in code ... in code first atom is N
 *  @todo CHECK!!
 */
#define CAT_ATnames {"H","N","CA","C","O","CB"}

/** @brief Number of atoms per residue */
#define NATOM 6

/** @brief Height of CB atoms relative to C-CA-N plane  */
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
 * sist->Rbond[5], used for Bending instead of the angle (for backward compatibility with published works).
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
 * @brief Main data structure for the caterpillar protein model
 *
 * This data structure stores the coordinates of all atoms in the protein, the torsion angles, 
 * number of contacts, and  type of each residue.
 */
typedef struct {
	size_t n_atom_per_res;  /**< Number of atoms per residue (6 if CB included otherwise 5). See the N_ATOM macro. */
	size_t n_atoms;         /**< Number of atoms in the protein. */
	size_t n_res;           /**< Number of residues in the protein. */
	double **coord;         /**< Array storing all atoms' coordinates. */
	double ** N;            /**< Slice of coord array storing the coordinates of N atoms. N[i] corresponds to the N atom in the i-th residue. */
	double **CA;            /**< Slice of coord array storing the coordinates of Ca atoms. Ca[i] corresponds to the Ca atom in the i-th residue. */
	double ** C;            /**< Slice of coord array storing the coordinates of C atoms. C[i] corresponds to the C atom in the i-th residue. */
	double ** O;            /**< Slice of coord array storing the coordinates of O atoms. O[i] corresponds to the O atom in the i-th residue. */
	double ** H;            /**< Slice of coord array storing the coordinates of H atoms. H[i] corresponds to the H atom in the i-th residue. */
	double **CB;            /**< Slice of coord array storing the coordinates of Cb atoms. Cb[i] corresponds to Cb atom in the i-th residue. */
	double *phi;            /**< Array storing phi torsion angles. phi[i] corresponds to the phi dihedral in the i-th residue. */
	double *psi;            /**< Array storing psi torsion angles. psi[i] corresponds to the psi dihedral in the i-th residue. */
	double *contacts;       /**< Array storing the number of contacts for each residue. This value is used to compute water exposure. */
	int    *residues;       /**< Array storing all residue types. Values are determined by #CAT_AACODE. Translation to one letter FASTA code is given in #CAT_FASTA. */
} cat_prot;

/**
 * @brief  Allocates and initializes a cat_prot structure
 *
 */
cat_prot * CAT_prot_alloc (size_t n_res, size_t n_atom_per_res);

/**
 * @brief Deallocates a cat_prot structure
 */
void CAT_prot_free (cat_prot * pr);

/**
 * @brief Prints atom coordinates on stdout
 */
void CAT_print(cat_prot *protein);

/**
 * @brief Constructs a protein configuration cat_prot from dihedral angles
 *
 */
cat_prot * CAT_build_from_dihed ( int n_res, size_t n_atom_per_res, double *orig, double *dihed, char *seq);

/**
 * @brief Constructs a liner protein configuration cat_prot
 *
 *
 */
void CAT_set_prot_linear ( cat_prot * protein, double *orig, double alpha );
/**
 * @brief Sets protein residue types from FASTA string
 */
void CAT_set_residues_fasta ( cat_prot * protein, int Seq_Length, char *Enc);

/**
 * @brief Adds hydrogens to an existing protein backbone
 */
void CAT_insert_hydrogens ( cat_prot * protein );

/**
 * @brief Adds a  CB atom to a given residue
 */
void CAT_insert_cbeta ( cat_prot *protein, int res_i, double azimut, double bond_length	);

/**
 * @brief Corrects numerical distortions in a cat_prot backbone
 */
void CAT_rescale ( cat_prot * protein );

/**
 * @brief Inserts peptide I based on the angles phi, psi, and alpha (link twist)
 *
 */
int CAT_add_peptide ( cat_prot *p, int I, double phi, double alpha, double psi );

/**
 * @brief Computes all dihedral angles
 *
 */
void CAT_prot_dihedrals (cat_prot *protein);

/**
 * @brief Computes the phi dihedral of a given residue c
 */
double compute_phi(cat_prot *p,int c);

/**
 * @brief Computes the psi dihedral of a given residue c
 */
double compute_psi(cat_prot *p,int c);

/**
 * @brief Recomputes all dihedral angles
 *
 * @todo have to be rewritten
 *
 * 
 */
void compute_dihedrals(cat_prot *p);

/**
 * @brief Calculates the torsion angle between four successive atoms on a backbone
 *
 */
double calc_dihedralf_angle(double *atom_1, double *atom_2, double *atom_3, double *atom_4);

/**
 * @brief Builds a caterpillar peptide in a local DH base.
 */
void build_peptide ( gsl_matrix *pep);
#endif //__CATERPILLAR__
