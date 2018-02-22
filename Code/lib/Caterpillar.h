#ifndef __CATERPILLAR__
#define __CATERPILLAR__

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_matrix.h>

//Alphabet
#define CAT_S 21
//Acceptable residues names. DUM for Dummy residue.
#define CAT_AADefs {"ALA","CYS", "ASP", "GLU", "PHE", "GLY", "HIS", "ILE",\
	"LYS" ,"LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL",\
	"TRP", "TYR","DUM"}
#define CAT_FASTA "ACDEFGHIKLMNPQRSTVWYX" 
//This differ from ivan's fasta decoder. The reason is that it makes it much 
//simpler to assign the code to the letter and viceversa (the num. id. is simply
//the index of the letter in the string).
#define CAT_AACODE {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20}
#define ATOM_H 0
#define ATOM_N 1
#define ATOM_CA 2
#define ATOM_C 3
#define ATOM_O 4
#define ATOM_CB 5
//Atom names
#define CAT_ATnames {"H","N","CA","C","O","CB"}
#define NATOM 6
#define CAT_Cb_AZIMUTH 1.2
#define CAT_Rbond_CCb 1.000000

//--Structure of residue
//Angles for H insertion (dot products actually)
#define CAT_Plat__C_HdotC_N (2.4318886938725948)		//angle HCN
#define CAT_Plat__N_CAdotC_N  (2.8923429607978419)		//angle CNCa
#define CAT_Plat__C_HdotC_O  (-2.1259814163093753)	//angle OCH
//Lengths
#define CAT_Rbond_CaN 1.4500000 //sist->Rbond[0]= 1.45000;
#define CAT_Rbond_CCa 1.5200000	//sist->Rbond[1]=1.52000;
#define CAT_Rbond_CO  1.2300000	//sist->Rbond[2]= 1.23000;
#define CAT_Rbond_CN  1.3300000	//sist->Rbond[3]= 1.33000;
#define CAT_Rbond_NH  1.0000000	//sist->Rbond[4]= 1.00000;
#define CAT_Rbond_NC  2.4479800	//sist->Rbond[5], used for Bending instead of the angle..
//Backbone angles
#define CAT_angle_NCaC 1.9373154697137058
#define CAT_angle_CaCN 2.017600615305445
#define CAT_angle_CNCa 2.1275563581810877
#define CAT_angle_CNH  2.0856684561332237
#define CAT_angle_CaCO 2.1135937241651326
//others



typedef struct {
	size_t n_atom_per_res;
	size_t n_atoms;
	size_t n_res;
	double **coord;
	double ** N;
	double **CA;
	double ** C;
	double ** O;
	double ** H;
	double **CB;
	double *phi;
	double *psi;
	double *contacts;
	int *residues;
} cat_prot;


cat_prot * CAT_prot_alloc (size_t n_res, size_t n_atom_per_res);
void CAT_prot_free (cat_prot * pr);
void CAT_print(cat_prot *protein);

cat_prot * CAT_build_linear_prot 	(	int n_res, size_t n_atom_per_res, double *orig, double alpha );
cat_prot * CAT_build_from_dihed 	( int n_res, size_t n_atom_per_res, double *orig, double *dihed, char *seq);
void CAT_set_prot_linear 					( cat_prot * protein, double *orig, double alpha );
void CAT_set_residues 						( cat_prot * protein, int *Dec);
void CAT_set_residues_fasta 			( cat_prot * protein, int Seq_Length, char *Enc);
void CAT_insert_hydrogens					( cat_prot * protein );
void CAT_insert_cbeta 						( cat_prot *protein, int res_i, double azimut, double bond_length	);
void CAT_rescale									( cat_prot * protein );
//void CAT_add_peptide 							( cat_prot * protein, int i, double phi, double psi);
//void CAT_add_peptide ( cat_prot * protein, int I, double phi, double psi, double angle_NCaC);

int 	CAT_add_peptide ( cat_prot *p, int I, double phi, double alpha, double psi );
//void CAT_add_peptide ( cat_prot * protein, int I, double phi, double psi, double angle_NCaC);
void CAT_prot_dihedrals (cat_prot *protein);
double compute_phi(cat_prot *p,int c);
double compute_psi(cat_prot *p,int c);

//these two have to be rewritten 
void compute_dihedrals(cat_prot *p);
double calc_dihedralf_angle(double *atom_1, double *atom_2, double *atom_3, double *atom_4);
#endif
