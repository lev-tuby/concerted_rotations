#include <stddef.h>
#include <stdlib.h>
#include "Caterpillar_energies.h"
#include "messages.h"
#include "my_geom.h"
#include "geom_prop.h"

// More general Saw function ...
inline double CATENR_Saw_g 	(double *r_atom_i, double *r_atom_j, double saw)
{
	return dist_d(r_atom_i, r_atom_j, 3) < saw ? CATENR_INFTY : CATENR_ZERO_EN;
}
inline double CATENR_Saw_eq_g 	(double r, double saw)
{
	return r<saw ? CATENR_INFTY : CATENR_ZERO_EN;
}
inline double CATENR_Saw 	(double *r_Ca_i, double *r_Ca_j)
{
	return dist_d(r_Ca_i,r_Ca_j,3)<CATENR_SAW_radius ? CATENR_INFTY : CATENR_ZERO_EN;
}
inline double CATENR_Saw_eq 	(double r)
{
	return r<CATENR_SAW_radius ? CATENR_INFTY : CATENR_ZERO_EN;
}

/*
double CATENR_Saw 	(double *r_Ca_i, double *r_Ca_j)
{
	if(dist_d(r_Ca_i,r_Ca_j,3)<CATENR_SAW_radius)
	{
		return CATENR_INFTY;
	}
	return CATENR_ZERO_EN;
}
*/

inline double CATENR_bond (double r)
{
	return 1.0/(1.0 + exp(2.5*(r-CATENR_CaCa_RANGE)));
}

inline double CATENR_Ca_Ca_eq	(int type_i, int type_j, double bond, double **M, double pref)
{
	return pref*bond*M[type_i][type_j];
}
double CATENR_Ca_Ca (double *r_Ca_i, int type_i, double *r_Ca_j, int type_j,int n, double **M, double pref)
{
	//LT tested OK 11/07/17 
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
double CATENR_HB  	(double *restrict r_H_i, double *restrict r_N_i, double *restrict r_O_j, double *restrict r_C_j, double pref)
{
	//LT Tested OK (?) 12/07/17. There is still a small discrepancy that I suppose is caused by pdb errors. must check with the binary file.
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
double CATENR_HB_eq  	(double *r_H_i, double *r_N_i, double *r_O_j, double *r_C_j,double r, double pref)
{
	//LT Tested OK (?) 12/07/17. There is still a small discrepancy that I suppose is caused by pdb errors. must check with the binary file.
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

inline double CATENR_Wat_eq		(int type, double contacts, double *Hydrophob,double *treshold,double pref)
{
	double ener;
	ener=Hydrophob[type]*(treshold[type]-contacts);
	return ener>0?ener:0;
}
double CATENR_Wat		(cat_prot *protein, int i, int n, double *Hydrophob, double *treshold, double pref)
{
	//LT Tested OK 11/07/17 ( reproduces energy of a protein in input. Mind the prefactors!)
	int j;
	double r;
	double bond;
	double ener;
	int type=protein->residues[i];
	bond=0;
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

double CATENR_Ca_Ca_Bend 	(double *r_Ca_i, double *r_Ca_j, double ref_bond, double pref)
{
	double r;
	r=dist_d(r_Ca_i,r_Ca_j,3);
	return pref*(r-ref_bond)*(r-ref_bond);
}


void CATENR_read_energy_params(int n, double **ptr_Hydrophob, double **ptr_H_treshold, double ***ptr_CaCa_Matrix, char filename[1024])
{
	//LT Tested OK 08/05/17 (read wpot file correctly)
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


