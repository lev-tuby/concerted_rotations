/**
 * @file
 * @brief File contain function for energy calculation
 * @todo It would be nice to take all energy functions from main compiled programs and move them here ...
 */

#include <stddef.h>
#include <stdlib.h>
#include "Caterpillar_energies.h"
#include "messages.h"
#include "my_geom.h"
#include "geom_prop.h"


/**
 * @brief Function calculate hard-core repulstion energy between two atoms
 *
 * @param[in]  *r_atom_i         Coordinates of first atom in format #d1t(3).
 * @param[in]  *r_atom_j         Coordinates of second atom in format #d1t(3).
 * @param[in]   saw              Minimal allowed distance between atoms.
 *
 * @return energy 0 or \f$ \infty \f$
 */
inline double CATENR_Saw (double *r_atom_i, double *r_atom_j, double saw)
{
    return dist_d(r_atom_i, r_atom_j, 3) < saw ? CATENR_INFTY : CATENR_ZERO_EN;
}

/**
 * @brief Function calculate hard-core repulstion energy from distance and cut-off
 *
 * @param[in]   r                Distance between atoms.
 * @param[in]   saw              Minimal allowed distance between atoms.
 *
 * @return energy 0 or \f$ \infty \f$
 */
inline double CATENR_Saw_eq 	(double r, double saw)
{
	return r<saw ? CATENR_INFTY : CATENR_ZERO_EN;
}

/**
 * @brief Function calculate attractive part of Catepillar interaction
 *
 * Function calculate continuous square well like potential.
 * By multipliing this value with Ca-Ca interaction energy we get whole Ca-Ca interaction.
 *
 * @todo I think name of the function is misleading ...
 * @todo Add shape of given potential in spplementary graph in documentation.
 *
 * @param[in]   r                Distance between atoms.
 *
 * @return non-scaled interaction energy
 */
inline double CATENR_bond (double r)
{
	return 1.0/(1.0 + exp(2.5*(r-CATENR_CaCa_RANGE)));
}

/**
 * @brief Function calculate attractive part of Catterpillar interaction
 *
 * Function take atom types of both Ca atoms, interaction matrix M and preffactor.
 * bond is calculated via #CATENR_bond.
 *
 * @param[in]   type_i           Residue type of first Ca atom.
 * @param[in]   type_j           Residue type of second Ca atom.
 * @param[in]   bond             Value of squre well like potentil based on Ca atom distance see #CATENR_bond.
 * @param[in] **M                Interaction matrix for all Ca-Ca interactions of different residues with dimension (#CAT_S-1) * (#CAT_S-1).
 * @param[in]   pref             Scaling prefactor of all Ca-Ca interactions #CATENR_CaCa_PREF.
 *
 * @return interaction energy
 */
inline double CATENR_Ca_Ca_eq	(int type_i, int type_j, double bond, double **M, double pref)
{
	return pref*bond*M[type_i][type_j];
}

/**
 * @brief All in one function for calculation of attractive part of Caterpillar model
 *
 * LT tested OK 11/07/17 
 * @note Is it right that prefactor pref is being used even if LINK_DIAGNOSE macro being invoked?
 * @todo remove n variable ... not used at all.
 *
 * @param[in]  *r_Ca_i           First Ca atom position in format d1t(3).
 * @param[in]   type_i           Residue type of first Ca atom.
 * @param[in]  *r_Ca_j           Second Ca atom position in format d1t(3).
 * @param[in]   type_j           Residue type of second Ca atom.
 * @param[in]   n                Well it should be size of interactio matrix **M ... but it is not used and size is defined by macro anyway
 * @param[in] **M                Interaction matrix for all Ca-Ca interactions of different residues with dimension (#CAT_S-1) * (#CAT_S-1).
 * @param[in]   pref             Scaling prefactor of all Ca-Ca interactions #CATENR_CaCa_PREF.
 *
 * @return interaction energy
 */
double CATENR_Ca_Ca (double *r_Ca_i, int type_i, double *r_Ca_j, int type_j,int n, double **M, double pref)
{
	double r;
	double ener=CATENR_ZERO_EN;
	r=dist_d(r_Ca_i,r_Ca_j,3);
#ifdef LINK_DIAGNOSE
	if(r<CATENR_CaCa_RANGE){
		ener=1;
	}else{
		ener=0;
	}
#else
	if(r<CATENR_CaCa_CUTOFF)
	{
		//ener=1.0/(1.0 + exp(2.5*(r-CATENR_CaCa_RANGE)));
		//ener*=M[type_i][type_j];
		ener=CATENR_bond(r)*M[type_i][type_j];
	}
#endif
	return pref*ener;
}

/*
double CATENR_HB  	(double *r_H_i, double *r_N_i, double *r_O_j, double *r_C_j,double pref)
{
	//LT tested OK 11/07/17 (reproduces energy of a protein given in input)
	// HB specific parameters
	double HB_dotprod_OHN_CUTOFF=0; //90 degrees angles
	double HB_dotprod_HOC_CUTOFF=0;
	double HB_power=2;
	double HB_LJ12=6;//EH_LJPower1_2 = 12/2
	double HB_LJ10=5;//EH_LJPower2_2 =10/2
	double HB_sigma12=pow(2.0,2*HB_LJ12)*HB_LJ10/(HB_LJ12-HB_LJ10);
	double HB_sigma10=pow(2.0,2*HB_LJ10)*HB_LJ12/(HB_LJ12-HB_LJ10);
	//-----
	int j;
	double r,v_HO[3];
	double r10,r12;
	double v_NH[3], v_OC[3];
	double dotprod_OHN;
	double dotprod_HOC;
	double ener=CATENR_ZERO_EN;

	r=dist_d(r_H_i,r_O_j,3);
	if(r<CATENR_HB_CUTOFF)
	{
		for(j=0;j<3;j++)
		{
			v_HO[j]=r_H_i[j]-r_O_j[j];
			v_NH[j]=r_N_i[j]-r_H_i[j];
			v_OC[j]=r_C_j[j]-r_O_j[j];
		}
		//normalize_d(&v_NH[0],3);
		//normalize_d(&v_OC[0],3);
		dotprod_OHN=-scal_d(&v_NH[0],&v_HO[0],3)/r;
		dotprod_HOC=scal_d(&v_OC[0],&v_HO[0],3)/(r*CAT_Rbond_CO);
		if(dotprod_OHN<=HB_dotprod_OHN_CUTOFF && dotprod_HOC<=HB_dotprod_HOC_CUTOFF)
		{
			r10 =pow(r,2*HB_LJ10 );
			r12=pow(r,2*HB_LJ12);
			ener=pow(dotprod_OHN*dotprod_HOC,HB_power);
			ener*=pref*(HB_sigma12/r12-HB_sigma10/r10);
			if (ener>100){ener=100;}
		}
	}
	return ener;
}
*/

/**
 * @brief Calculate H-bond interaction from positions of N H and O atoms
 *
 * LT Tested OK (?) 12/07/17. There is still a small discrepancy that I suppose is caused by pdb errors. must check with the binary file.
 * @todo Since prefactor is given as macro it migh be removed from parameters ... unless same function being used for optimization purposes...
 *
 * @param[in]  *r_H_i            First Ca atom position in format d1t(3).
 * @param[in]  *r_N_i            Residue type of first Ca atom.
 * @param[in]  *r_O_j            Second Ca atom position in format d1t(3).
 * @param[in]  *r_C_j            Residue type of second Ca atom.
 * @param[in]   pref             Scaling prefactor of all HB interactions #CATENR_HB_PREF.
 *
 * @return HB interaction energy
 */
double CATENR_HB (double *restrict r_H_i, double *restrict r_N_i, double *restrict r_O_j, double *restrict r_C_j, double pref)
{
	// HB specific parameters
	double HB_dotprod_OHN_CUTOFF=0; //90 degrees angles
	double HB_dotprod_HOC_CUTOFF=0;
	double HB_power=2;
	//double HB_LJ12=6;//EH_LJPower1_2 = 12/2
	//double HB_LJ10=5;//EH_LJPower2_2 =10/2
	//double HB_sigma12=pow(2.0,2*HB_LJ12)*HB_LJ10/(HB_LJ12-HB_LJ10);
	//double HB_sigma10=pow(2.0,2*HB_LJ10)*HB_LJ12/(HB_LJ12-HB_LJ10);
	double HB_LJ12=12;//EH_LJPower1_2 = 12/2
	double HB_LJ10=10;//EH_LJPower2_2 =10/2
	double HB_LJ12_pref=5;//EH_LJPower1_2 = 12/2
	double HB_LJ10_pref=6;//EH_LJPower2_2 =10/2
	double HB_sigma12=pow(2.0,HB_LJ12);
	double HB_sigma10=pow(2.0,HB_LJ10);
	//-----
	int j;
	double r,v_HO[3];
	double r10,r12;
	double v_NH[3], v_OC[3];
	double dotprod_OHN;
	double dotprod_HOC;
	double ener=CATENR_ZERO_EN;

	static int nbond=0;
	r=dist_d(r_H_i,r_O_j,3);
	if(r<CATENR_HB_CUTOFF)
	{
		for(j=0;j<3;j++)
		{
			v_HO[j]=r_H_i[j]-r_O_j[j];
			v_NH[j]=r_N_i[j]-r_H_i[j];
			v_OC[j]=r_C_j[j]-r_O_j[j];
		}
		normalize_d(&v_NH[0],3);
		normalize_d(&v_OC[0],3);
		dotprod_OHN=-scal_d(&v_NH[0],&v_HO[0],3)/r;
		dotprod_HOC=scal_d(&v_OC[0],&v_HO[0],3)/r;
		//fprintf(stderr,"OHN: %10.5f COH: %10.5f r %10.5f\n",acos(dotprod_OHN)*180.0/M_PI,acos(dotprod_HOC)*180.0/M_PI,r);
		if(dotprod_OHN<=HB_dotprod_OHN_CUTOFF && dotprod_HOC<=HB_dotprod_HOC_CUTOFF)
		{
			r10 =pow(r,HB_LJ10);
			r12=pow(r,HB_LJ12);
			ener=pow(dotprod_OHN*dotprod_HOC,HB_power);
			ener*=(HB_LJ12_pref*HB_sigma12/r12-HB_LJ10_pref*HB_sigma10/r10);
			ener*=pref;
			if (ener>100){ener=100;}
		}
	}
	return ener;
}

/**
 * @brief Calculate H-bond interaction from positions of N H and O atoms
 *
 * LT Tested OK (?) 12/07/17. There is still a small discrepancy that I suppose is caused by pdb errors. must check with the binary file.
 * @note This function in comparison to #CATENR_HB do not check if atoms are within HB cut-off distance #CATENR_HB_CUTOFF.
 * @todo Since prefactor is given as macro it migh be removed from parameters ... unless same function being used for optimization purposes...
 *
 * @param[in]  *r_H_i            First Ca atom position in format d1t(3).
 * @param[in]  *r_N_i            Residue type of first Ca atom.
 * @param[in]  *r_O_j            Second Ca atom position in format d1t(3).
 * @param[in]  *r_C_j            Residue type of second Ca atom.
 * @param[in]   pref             Scaling prefactor of all HB interactions #CATENR_HB_PREF.
 *
 * @return HB interaction energy
 */
double CATENR_HB_eq  	(double *r_H_i, double *r_N_i, double *r_O_j, double *r_C_j,double r, double pref)
{
	// HB specific parameters
	double HB_dotprod_OHN_CUTOFF=0; //90 degrees angles
	double HB_dotprod_HOC_CUTOFF=0;
	double HB_power=2;
	double HB_LJ12=12;//EH_LJPower1_2 = 12/2
	double HB_LJ10=10;//EH_LJPower2_2 =10/2
	double HB_LJ12_pref=5;//EH_LJPower1_2 = 12/2
	double HB_LJ10_pref=6;//EH_LJPower2_2 =10/2
	double HB_sigma12=pow(2.0,HB_LJ12);
	double HB_sigma10=pow(2.0,HB_LJ10);
	//-----
	int j;
	double v_HO[3];
	double r10,r12;
	double v_NH[3], v_OC[3];
	double dotprod_OHN;
	double dotprod_HOC;
	double ener=CATENR_ZERO_EN;

	static int nbond=0;
	//r=dist_d(r_H_i,r_O_j,3);
	for(j=0;j<3;j++)
	{
		v_HO[j]=r_H_i[j]-r_O_j[j];
		v_NH[j]=r_N_i[j]-r_H_i[j];
		v_OC[j]=r_C_j[j]-r_O_j[j];
	}
	normalize_d(&v_NH[0],3);
	normalize_d(&v_OC[0],3);
	dotprod_OHN=-scal_d(&v_NH[0],&v_HO[0],3)/r;
	dotprod_HOC=scal_d(&v_OC[0],&v_HO[0],3)/r;
	//fprintf(stderr,"OHN: %10.5f COH: %10.5f r %10.5f\n",acos(dotprod_OHN)*180.0/M_PI,acos(dotprod_HOC)*180.0/M_PI,r);
	if(dotprod_OHN<=HB_dotprod_OHN_CUTOFF && dotprod_HOC<=HB_dotprod_HOC_CUTOFF)
	{
		r10 =pow(r,HB_LJ10);
		r12=pow(r,HB_LJ12);
		ener=pow(dotprod_OHN*dotprod_HOC,HB_power);
		ener*=(HB_LJ12_pref*HB_sigma12/r12-HB_LJ10_pref*HB_sigma10/r10);
		ener*=pref;
		if (ener>100){ener=100;}
	}
	return ener;
}

/**
 * @brief Function calculate Caterpillar implicit solvatation energy
 *
 * Function substract number of contacts of given residue with other residues from treshold.
 * Treshold value being larger for hydrophobic residues, which if not burried (large number of contacts) would produce
 * energy penalty scaled by hydrophopbicity of residue *Hydrophob.
 *
 * @param[in]   type             Residue type one of #CAT_AACODE.
 * @param[in]   contacts         Number of residues within interaction range (corespond to degree of residue burrial).
 * @param[in]  *Hydrophob        Hydrophobicity energy scaling.
 * @param[in]  *treshold         Number of residue contacts that have to be satisfied otherwise energy penalty is produced.
 * @param[in]   pref             Scaling prefactor for all solvatation energy.
 *
 * @return hydrophobic interaction energy
 */
inline double CATENR_Wat_eq (int type, double contacts, double *Hydrophob,double *treshold,double pref)
{
	double ener;
	ener=Hydrophob[type]*(treshold[type]-contacts);
	return ener>0?ener:0;
}

/**
 * @brief Function calculate Caterpillar implicit solvatation energy of one residue
 *
 * Function substract number of contacts of given residue bond with other residues from treshold.
 * Treshold value being larger for hydrophobic residues, which if not burried (large number of contacts) would produce
 * energy penalty scaled by hydrophopbicity of residue *Hydrophob.
 * LT Tested OK 11/07/17 ( reproduces energy of a protein in input. Mind the prefactors!)
 *
 * @note Function iterate only on given protein chain, so that in case of multiple proteins
 * this will not calculate solvation energy of whole system but only within one protein!
 *
 * @todo Remove not used value n.
 *
 * @param[in]  *protein          Protein backbone representation.
 * @param[in]   i                Index of residue to calculate its solvatation energy.
 * @param[in]   n                Not used value.
 * @param[in]  *Hydrophob        Hydrophobicity energy scaling.
 * @param[in]  *treshold         Number of residue contacts that have to be satisfied otherwise energy penalty is produced.
 * @param[in]   pref             Scaling prefactor for all solvatation energy.
 *
 * @return hydrophobic energy of residue i
 */
double CATENR_Wat (cat_prot *protein, int i, int n, double *Hydrophob, double *treshold, double pref)
{
	int j;
	double r;
	double bond;
	double ener;
	int type=protein->residues[i];
	bond=0; // number of all contacts of residue i
	for(j=0;j<i-2;j++)
	{
		r=dist_d(protein->CA[i],protein->CA[j],3);
		if(r<CATENR_CaCa_CUTOFF)
		{
			//bond+=1.0/(1.0 + exp(2.5*(r-CATENR_CaCa_RANGE)));
			bond+=CATENR_bond(r);
		}
	}
	for(j=i+3;j<protein->n_res;j++)
	{
		r=dist_d(protein->CA[i],protein->CA[j],3);
		if(r<CATENR_CaCa_CUTOFF)
		{
			//bond+=1.0/(1.0 + exp(2.5*(r-CATENR_CaCa_RANGE)));
			bond+=CATENR_bond(r);
		}
	}
	ener=Hydrophob[type]*(treshold[type]-bond);
	ener=ener>0?ener:0;
	protein->contacts[i]=bond;
	return pref*ener;
}

/**
 * @brief Function calculate Caterpillar implicit solvatation energy of all residue
 *
 * @note Function iterate only on given protein chain, so that in case of multiple proteins
 * this will not calculate solvation energy of whole system but only within one protein!
 *
 * @param[in]  *protein          Protein backbone representation.
 * @param[in]  *Hydrophob        Hydrophobicity energy scaling.
 * @param[in]  *treshold         Number of residue contacts that have to be satisfied otherwise energy penalty is produced.
 * @param[in]   pref             Scaling prefactor for all solvatation energy.
 *
 * @return hydrophobic energy of protein
 */
double CATENR_Wat_total		(cat_prot *protein, double *Hydrophob, double *treshold, double pref)
{
	int i,j;
	double r;
	double bond;
	double ener;
	double Ewat=0;
	int type;
	for(i=0;i<protein->n_res;i++)
	{
		bond=0;
		type=protein->residues[i];
		for(j=0;j<i-2;j++)
		{
			r=dist_d(protein->CA[i],protein->CA[j],3);
			if(r<CATENR_CaCa_CUTOFF)
			{
				//bond+=1.0/(1.0 + exp(2.5*(r-CATENR_CaCa_RANGE)));
				bond+=CATENR_bond(r);
			}
		}
		for(j=i+3;j<protein->n_res;j++)
		{
			r=dist_d(protein->CA[i],protein->CA[j],3);
			if(r<CATENR_CaCa_CUTOFF)
			{
				//bond+=1.0/(1.0 + exp(2.5*(r-CATENR_CaCa_RANGE)));
				bond+=CATENR_bond(r);
			}
		}
		ener=Hydrophob[type]*(treshold[type]-bond);
		ener=ener>0?ener:0;
		protein->contacts[i]=bond;
		Ewat+=pref*ener;
	}
	return Ewat;
}

/**
 * @brief Function calculate bending energy
 *
 *
 * @param[in]  *r_Ca_i           Position of first Ca atom #d1t(3).
 * @param[in]  *r_Ca_j           Position of second Ca atom #d1t(3).
 * @param[in]   ref_bond         Reference length between two Ca atoms.
 * @param[in]   pref             Scaling prefactor for bending energy.
 *
 * @return hydrophobic energy of protein
 */
double CATENR_Ca_Ca_Bend 	(double *r_Ca_i, double *r_Ca_j, double ref_bond, double pref)
{
	double r;
	r=dist_d(r_Ca_i,r_Ca_j,3);
	return pref*(r-ref_bond)*(r-ref_bond);
}

/**
 * @brief Function read energy parameters from text file
 *
 * Function read Ca-Ca attractive interaction, hydrophobic scale, and hydrophobic treshold.
 *
 * LT Tested OK 08/05/17 (read wpot file correctly)
 *
 * @param[in]        n                Number of different residue types (#CAT_S-1).
 * @param[in,out]  **ptr_Hydrophob    Pointer to hydrophobic scale in format d1t(n).
 * @param[in,out]  **ptr_H_treshold   Pointer to hydrophobic treshold in format d1t(n).
 * @param[in,out] ***ptr_CaCa_Matrix  Symmetric matrix of Ca-Ca interaction between different residue types in format #d2t(n,n).
 * @param[in]        filename         File with potential.
 *
 * @return hydrophobic energy of protein
 */
void CATENR_read_energy_params(int n, double **ptr_Hydrophob, double **ptr_H_treshold, double ***ptr_CaCa_Matrix, char filename[1024])
{
	double *Hydrophob	=d1t(n);
	double *H_treshold =d1t(n);
	double **CaCa_Matrix=d2t(n,n);
	FILE   *fp;
	fp=fopen(filename,"r");
	if ( fp == NULL) { failed("File not found!\n"); }
	rdvector		(n,	H_treshold,	fp);
	rdvector		(n,	Hydrophob,	fp);
	rdarray_symm(n,	CaCa_Matrix,fp);

	*ptr_H_treshold=H_treshold;
	*ptr_Hydrophob=Hydrophob;
	*ptr_CaCa_Matrix=CaCa_Matrix;
	fclose(fp);
}


