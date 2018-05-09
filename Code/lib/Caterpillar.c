/**
 * @file
 * @brief Function for manipulation of Cterpillar peptide backbone
 * @todo There are some large comented regions which might be good to remove if not relevant anymore.
 * @todo Some functions might still need refinement ... regarding how they are written.
 * @todo CAT_prot_dihedrals and compute_dihedrals does same thing ... I would leave only one ...
 * @todo ... too many functions for dihedral calculation ... it would be best to have one that calcualte from 4 atoms then one that does from rezidue in protein and one that calculate for all residues
 * @note What about gly residue and CB atom ... how to treat that ... keep CB there but do not use it probbably ...
 */

#include "Caterpillar.h"
#include "./messages.h"
#include "./my_geom.h"
#include "./geom_prop.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include "./quaternions.h"

/**
 * @brief prototypes for internal use
 *
 */
void build_peptide ( gsl_matrix *pep);

/**
 * @brief Function initialize cat_prot strucure
 *
 * All coordinates are set to 0.0, all dihedrals are set to 0.0. p->=n_res*n_atom_per_res, p->n_res=n_res, p->n_atom_per_res=n_atom_per_res, p->residues[i]=-1, p->contacts[i]=-1.0.
 * Atoms are aliased ... it is kind of dangerous though. It remains unclear how optimization is affected inside a function.
 * 
 * @param[in]        n_res           Number of residues in peptide
 * @param[in]        n_atom_per_res  Number of atoms in each residue
 *
 * @return initialized empty *cat_prot
 */
cat_prot * CAT_prot_alloc (size_t n_res, size_t n_atom_per_res)
{
    int
        i;
    double
        *a = NULL;

    size_t
        n;
    cat_prot
        *p = NULL;

    if ((p = (cat_prot *)malloc(sizeof(cat_prot))) == NULL)
    {
        failed("CAT_prot_alloc failed allocating p");
    }
    //Allocate  arrays
    n = n_res*n_atom_per_res;
    if ((p->residues  = (int *) malloc (n_res * sizeof (int ))) == NULL)
        failed ("CAT_prot_alloc: failed residues");
    if ((p->contacts  = (double *) calloc (n_res, sizeof (double))) == NULL)
        failed ("CAT_prot_alloc: failed contacts");
    if ((p->phi  = (double *) malloc (n_res * sizeof (double))) == NULL)
        failed ("CAT_prot_alloc: failed phi");
    if ((p->psi  = (double *) malloc (n_res * sizeof (double))) == NULL)
        failed ("CAT_prot_alloc: failed psi");
    if ((p->coord = (double **) malloc (n * sizeof (double *))) == NULL)
        failed ("CAT_prot_alloc: failed coord");
    //pointer aliases
    if ((p->N  = (double **) malloc (n_res * sizeof (double *))) == NULL)
        failed ("CAT_prot_alloc: failed N");
    if ((p->CA = (double **) malloc (n_res * sizeof (double *))) == NULL)
        failed ("CAT_prot_alloc: failed CA");
    if ((p->C  = (double **) malloc (n_res * sizeof (double *))) == NULL)
        failed ("CAT_prot_alloc: failed C");
    if ((p->O  = (double **) malloc (n_res * sizeof (double *))) == NULL)
        failed ("CAT_prot_alloc: failed O");
    if ((p->H  = (double **) malloc (n_res * sizeof (double *))) == NULL)
        failed ("CAT_prot_alloc: failed H");
    if (n_atom_per_res==6)
    {
        if ((p->CB = (double **) malloc (n_res * sizeof (double *))) == NULL)
        failed ("CAT_prot_alloc: failed CB");
    }
    else
    {
        p->CB = NULL;
    }
    //set up a 2D array
    if ((p->coord[0] = (double *) malloc (n*3 * sizeof (double))) == NULL)
        failed ("d2t: failed n*3");
    for (i = 0; i < n - 1; i++)
        p->coord[i + 1] = p->coord[i] + 3;
    //set all coords to zero
    for (i = 0, a = p->coord[0]; i < n *3; i++)
        *a++ = 0;
    //alias the atoms...kind of dangerous though.
    //It remains unclear how optimization is affected inside a function.
    p->N[0] =p->coord[ATOM_N];
    p->CA[0]=p->coord[ATOM_CA];
    p->C[0] =p->coord[ATOM_C];
    p->O[0] =p->coord[ATOM_O];
    p->H[0] =p->coord[ATOM_H];
    if(n_atom_per_res==6)
    {
        p->CB[0]=p->coord[ATOM_CB];
    }
    for(i=0; i<n_res-1; i++)
    {
        p->CA[i+1] = p->CA[i] + 3 * n_atom_per_res;
        p->C [i+1] = p->C [i] + 3 * n_atom_per_res;
        p->O [i+1] = p->O [i] + 3 * n_atom_per_res;
        p->H [i+1] = p->H [i] + 3 * n_atom_per_res;
        p->N [i+1] = p->N [i] + 3 * n_atom_per_res;
        if(n_atom_per_res == 6)
        {
            p->CB[i+1] = p->CB[i] + 3 * n_atom_per_res;
        }
    }

    p->n_atoms        = n;
    p->n_res          = n_res;
    p->n_atom_per_res = n_atom_per_res;
    //set all residues and contacts to -1, all angles  to zero
    for(i = 0; i < n_res; i++)
    {
        p->residues[i] = -1  ;
        p->contacts[i] = -1.0;
        p->phi[i]      =  0.0;
        p->psi[i]      =  0.0;
    }
    return p;
}

/**
 * @brief Deallocation of cat_prot structure
 *
 * @param[in,out]    *p  cat_prot structure to dealloc
 *
 * @return \c void
 */
void CAT_prot_free ( cat_prot * p)
{
	if(p!=NULL)
	{
		free(p->coord[0]);
		free(p->coord);
		free(p->N);
		free(p->CA);
		free(p->C);
		free(p->O);
		free(p->H);
		if(p->CB!=NULL)
			free(p->CB);
		free(p->residues);
		free(p->contacts);
		free(p->phi);
		free(p->psi);
		free(p);
	}
}

/**
 * @brief Function print in stdout position of all atoms in protein
 *
 * 
 * @param[in]      *protein        Protein configuration is printed out
 *
 * @return \c void
 */
void CAT_print(cat_prot *protein)
{
    for(int i=0; i<protein->n_res;i++)
    {
        printf("Resi: %i\n", i);
        printf("N:  %g %g %g\n", protein->N[i][0], protein->N[i][1], protein->N[i][2]);
        printf("CA: %g %g %g\n", protein->CA[i][0], protein->CA[i][1], protein->CA[i][2]);
        printf("C:  %g %g %g\n", protein->C[i][0], protein->C[i][1], protein->C[i][2]);
        printf("O:  %g %g %g\n", protein->O[i][0], protein->O[i][1], protein->O[i][2]);
        printf("H:  %g %g %g\n", protein->H[i][0], protein->H[i][1], protein->H[i][2]);
        if(protein->n_atom_per_res == 6)
        {
            printf("CB: %g %g %g\n", protein->CB[i][0], protein->CB[i][1], protein->CB[i][2]);
        }
    }
}

/**
 * @brief Function build up cat_prot from dihedral angles
 *
 * Build a Caterpillar protein, assuming that the first dihedral angle passed through dihed is a phi angle.
 * Function is bit dangerous ... using n_atom_per_res other then 5 or 6 cause silent errors ... check it maybe.
 * 
 * @param[in]        n_res          Number of residues in protein
 * @param[in]        n_atom_per_res Number of atoms in one residue
 * @param[in]       *orig           Position of first atom (so far it is N)
 * @param[in]       *dihed          Array of protein dihedral angles (in radians) starting with phi
 * @param[in]       *seq            Protein sequnce in one letter code
 *
 * @return builded cat_prot
 */
cat_prot * CAT_build_from_dihed ( int n_res, size_t n_atom_per_res, double *orig, double *dihed, char *seq)
{
	cat_prot *p = CAT_prot_alloc( n_res, n_atom_per_res );
	CAT_set_prot_linear(p,orig,0.0);
	CAT_set_residues_fasta(p,n_res, seq);
	for(int i = 0; i < n_res; i++)
	{
		CAT_add_peptide(p,i,dihed[2*i],M_PI-CAT_angle_NCaC,dihed[2*i+1]);
	}
	return p;
}

/**
 * @brief Function build up liner protein chain
 *
 * Tested OK. LT 02.08.16
 * 
 * @param[in,out]   *protein        cat_prot structure to be remodeled
 * @param[in]        orig[3]        Starting poinf from where peptide chain is build
 * @param[in]        alpha          Deviation of 
 *
 * @return \c void
 */
void CAT_set_prot_linear ( cat_prot *protein, double orig[3], double alpha )
{
	size_t n_res = protein->n_res;
	//size_t n_atoms = protein->n_atoms;
	//size_t n_atom_per_res = protein->n_atom_per_res;
	size_t i;
	double theta;

	protein->N[0][0]=orig[0];
	protein->N[0][1]=orig[1];
	protein->N[0][2]=orig[2];
	//proceed two residues at a time, to set them in trans configurations
	for(i=0;i<n_res;i+=2)
	{
		//set CA i
		theta=((90.0-alpha)/180.0)*M_PI;
		protein->CA[i][0]=protein->N[i][0]+CAT_Rbond_CaN*cos(theta);
		protein->CA[i][1]=protein->N[i][1]+CAT_Rbond_CaN*sin(theta);
		protein->CA[i][2]=protein->N[i][2];
		//set C i
		theta=((270.0+111.0-alpha)/180.0)*M_PI;
		//theta+=(111.0-alpha)/180.0*M_PI-M_PI;
		protein->C[i][0]=protein->CA[i][0]+CAT_Rbond_CCa*cos(theta);
		protein->C[i][1]=protein->CA[i][1]+CAT_Rbond_CCa*sin(theta);
		protein->C[i][2]=protein->CA[i][2];
		//set O i
		//theta+=(111.0-alpha)/180.0*M_PI-M_PI;
		theta=((90.0+111.0+121.1-alpha)/180.0)*M_PI;
		protein->O[i][0]=protein->C[i][0]+CAT_Rbond_CO*cos(theta);
		protein->O[i][1]=protein->C[i][1]+CAT_Rbond_CO*sin(theta);
		protein->O[i][2]=protein->C[i][2];
		//set H i
		theta=((90.0+121.9+119.5-alpha)/180.0)*M_PI;
		protein->H[i][0]=protein->N[i][0]+CAT_Rbond_NH*cos(theta);
		protein->H[i][1]=protein->N[i][1]+CAT_Rbond_NH*sin(theta);
		protein->H[i][2]=protein->N[i][2];
		if(i<n_res-1)
		{
			//set N i+1
			theta=((90.0-115.6+111.0-alpha)/180.0)*M_PI;
			protein->N[i+1][0]=protein->C[i][0]+CAT_Rbond_CN*cos(theta);
			protein->N[i+1][1]=protein->C[i][1]+CAT_Rbond_CN*sin(theta);
			protein->N[i+1][2]=protein->C[i][2];
			//set CA i+1
			theta=((270.0+121.9-115.6+111.0-alpha)/180.0)*M_PI;
			protein->CA[i+1][0]=protein->N[i+1][0]+CAT_Rbond_CaN*cos(theta);
			protein->CA[i+1][1]=protein->N[i+1][1]+CAT_Rbond_CaN*sin(theta);
			protein->CA[i+1][2]=protein->N[i+1][2];
			//set C i+1
			theta=((90.0+(121.9-115.6)-alpha)/180.0)*M_PI;
			protein->C[i+1][0]=protein->CA[i+1][0]+CAT_Rbond_CCa*cos(theta);
			protein->C[i+1][1]=protein->CA[i+1][1]+CAT_Rbond_CCa*sin(theta);
			protein->C[i+1][2]=protein->CA[i+1][2];
			//set O i+1
			theta=((123.2-90.0+121.9-alpha)/180.0)*M_PI;
			protein->O[i+1][0]=protein->C[i+1][0]+CAT_Rbond_CO*cos(theta);
			protein->O[i+1][1]=protein->C[i+1][1]+CAT_Rbond_CO*sin(theta);
			protein->O[i+1][2]=protein->C[i+1][2];
			//set H i+1
			theta=((-90.0+118.2+(121.9-115.6+111.0)-alpha)/180.0)*M_PI;
			protein->H[i+1][0]=protein->N[i+1][0]+CAT_Rbond_NH*cos(theta);
			protein->H[i+1][1]=protein->N[i+1][1]+CAT_Rbond_NH*sin(theta);
			protein->H[i+1][2]=protein->N[i+1][2];
		}
		//set N i+2
		if(i<n_res-2)
		{
			theta=((270.0+121.9-alpha)/180.0)*M_PI;
			protein->N[i+2][0]=protein->C[i+1][0]+CAT_Rbond_CN*cos(theta);
			protein->N[i+2][1]=protein->C[i+1][1]+CAT_Rbond_CN*sin(theta);
			protein->N[i+2][2]=protein->C[i+1][2];
		}
	}
	for(i=0;i<n_res;i++)
	{
		protein->phi[i]=M_PI;
		protein->psi[i]=M_PI;
	}
	if(protein->CB!=NULL)
	{
		for(i=0;i<n_res;i++)
		{
			CAT_insert_cbeta(protein,i,CAT_Cb_AZIMUTH, CAT_Rbond_CCb);
		}
	}
}

/**
 * @brief Function set protein residue types from FASTA string
 *
 * 
 * @param[in,out]   *protein        Protein modified to corespond to FASTA sequnce
 * @param[in]        Seq_Length     Number of residues in FASTA formated string (*Enc)
 * @param[in]       *Enc            String of chars coresponding to FASTA code
 *
 * @return \c void
 */
void CAT_set_residues_fasta ( cat_prot * protein, int Seq_Length, char *Enc)
{
	int i,j;
	char fasta_code[CAT_S]=CAT_FASTA;
	int cat_aacode[CAT_S]=CAT_AACODE;


	if(Seq_Length != protein->n_res)
	{
		failed("CAT_set_residues_fasta: Seq_Length and n_res differ.");
	}
	for(i=0;i<Seq_Length;i++)
	{
		for(j=0;j<CAT_S;j++)
		{
			if(Enc[i]==fasta_code[j]) { protein->residues[i]=cat_aacode[j]; }
		}
	}
}

/**
 * @brief Function add hydrogens to existing peptide backbone
 *
 * 
 * @param[in,out]   *protein        cat_prot structure to be remodeled
 *
 * @return \c void
 */
void CAT_insert_hydrogens ( cat_prot * protein )
{
	int i;
	double a[3],b[3],c[3],prod[3],temp[3];
	double M[9];
	double theta;
	char error[2048];

	theta=((90.0+121.9+119.5-45.0)/180.0)*M_PI;
	protein->H[0][0]=protein->N[0][0]+CAT_Rbond_NH*cos(theta);
	protein->H[0][1]=protein->N[0][1]+CAT_Rbond_NH*sin(theta);
	protein->H[0][2]=protein->N[0][2];

	for(i=1;i<protein->n_res;i++)
	{
		a[0]=protein->O[i-1][0]-protein->C[i-1][0];
		a[1]=protein->O[i-1][1]-protein->C[i-1][1];
		a[2]=protein->O[i-1][2]-protein->C[i-1][2];

		b[0]=protein->N[i][0]-protein->C[i-1][0];
		b[1]=protein->N[i][1]-protein->C[i-1][1];
		b[2]=protein->N[i][2]-protein->C[i-1][2];

		// Cross prod
		prod[0]=(a[1]*b[2]-a[2]*b[1]);
		prod[1]=(a[2]*b[0]-a[0]*b[2]);
		prod[2]=(a[0]*b[1]-a[1]*b[0]);

		//Check conditions

		temp[0]=0; //Planarity
		temp[1]=CAT_Plat__C_HdotC_N;
		temp[2]=CAT_Plat__C_HdotC_O;

		M[0]=prod[0];
		M[1]=prod[1];
		M[2]=prod[2];
		M[3]=b[0];
		M[4]=b[1];
		M[5]=b[2];
		M[6]=a[0];
		M[7]=a[1];
		M[8]=a[2];
		//printf("iter. %d\n",i);
		//for(j=0;j<3;j++){
		//	printf(" %lf %lf %lf\n",M[j*3],M[j*3+1],M[j*3+2]);
		//}
		gsl_matrix_view A = gsl_matrix_view_array (M, 3, 3);
		gsl_vector_view v_b = gsl_vector_view_array (temp, 3);
 		gsl_vector *x = gsl_vector_alloc (3);
		gsl_linalg_HH_solve ( &A.matrix,  &v_b.vector, x);

		c[0]=gsl_vector_get (x,0);
		c[1]=gsl_vector_get (x,1);
		c[2]=gsl_vector_get (x,2);

		if(fabs(prod[0]*c[0]+prod[1]*c[1]+prod[2]*c[2]-temp[0])>1e-10){
			sprintf(error,"CAT_insert_Hs. Hydrogen %d. Violates planarity condition.\n",i);
			sprintf(error,"1 %10.5lf \n2 %10.5lf \n3 %10.5lf \n",prod[0]*c[0]+prod[1]*c[1]+prod[2]*c[2],b[0]*c[0]+b[1]*c[1]+b[2]*c[2],a[0]*c[0]+a[1]*c[1]+a[2]*c[2]);
			failed(error);
		}
		if(fabs(b[0]*c[0]+b[1]*c[1]+b[2]*c[2]-temp[1])>1e-10){
			sprintf(error,"CAT_insert_Hs. Hydrogen %d. Violates HCN angle.\n",i);
			sprintf(error,"1 %10.5lf \n2 %10.5lf \n3 %10.5lf \n",prod[0]*c[0]+prod[1]*c[1]+prod[2]*c[2],b[0]*c[0]+b[1]*c[1]+b[2]*c[2],a[0]*c[0]+a[1]*c[1]+a[2]*c[2]);
			failed(error);
		}
		if(fabs(a[0]*c[0]+a[1]*c[1]+a[2]*c[2]-temp[2])>1e-10){
			sprintf(error,"CAT_insert_Hs. Hydrogen %d. Violates OCH angle.\n",i);
			sprintf(error,"1 %10.5lf \n2 %10.5lf \n3 %10.5lf \n",prod[0]*c[0]+prod[1]*c[1]+prod[2]*c[2],b[0]*c[0]+b[1]*c[1]+b[2]*c[2],a[0]*c[0]+a[1]*c[1]+a[2]*c[2]);
			failed(error);
		}

		protein->H[i][0]=c[0]+protein->C[i-1][0];
		protein->H[i][1]=c[1]+protein->C[i-1][1];
		protein->H[i][2]=c[2]+protein->C[i-1][2];

		gsl_vector_free(x);
	}
}

/**
 * @brief Function add pozition o CB atom to given residue of protein
 *
 * @todo Check CaCb_versor calculation ... it seems that in the end it is kust -azimut*normal ... and mid part is unimportatnt? since CaCb_versor=-CaCb_versor-azimut*normal???
 *
 * @param[in,out]   *protein        Protein in which CB are added
 * @param[in]        res_i          Residue number to which CB is added
 * @param[in]        azimut         Azimut angle of CB bond
 * @param[in]        bond_length    Length of CA-CB bond
 *
 * @return \c void
 */
void CAT_insert_cbeta 	( cat_prot *protein, int res_i, double azimut, double bond_length	)
{
	double CaC[3];
	double CaN[3];
	double normal[3];
	double CaCb_versor[3];
	cat_prot *p=protein;

	CaC[0]=p->C[res_i][0]-p->CA[res_i][0];
	CaC[1]=p->C[res_i][1]-p->CA[res_i][1];
	CaC[2]=p->C[res_i][2]-p->CA[res_i][2];

	CaN[0]=p->N[res_i][0]-p->CA[res_i][0];
	CaN[1]=p->N[res_i][1]-p->CA[res_i][1];
	CaN[2]=p->N[res_i][2]-p->CA[res_i][2];

	vecprod_d(CaC,CaN,normal);
	normalize_d(normal,3);


	CaCb_versor[0]=CaC[0]+CaN[0];
	CaCb_versor[1]=CaC[1]+CaN[1];
	CaCb_versor[2]=CaC[2]+CaN[2];
	normalize_d(CaCb_versor,3);

	//This is the versor from the CA to the CB
	CaCb_versor[0]=-CaCb_versor[0]-azimut*normal[0];
	CaCb_versor[1]=-CaCb_versor[1]-azimut*normal[1];
	CaCb_versor[2]=-CaCb_versor[2]-azimut*normal[2];
	normalize_d(CaCb_versor,3);

	p->CB[res_i][0]=p->CA[res_i][0]+bond_length*CaCb_versor[0];
	p->CB[res_i][1]=p->CA[res_i][1]+bond_length*CaCb_versor[1];
	p->CB[res_i][2]=p->CA[res_i][2]+bond_length*CaCb_versor[2];
}
/*
void CAT_rescale( cat_prot *protein)
{
	int i,j;
	double b,bond[3];
	for(j=0;j<3;j++) { bond[j]=protein->H[0][j]-protein->N[0][j]; }
	b=norm_d(&(bond[0]),3);
	if(fabs(b-CAT_Rbond_NH)>1e-12)
	{
		for(j=0;j<3;j++)
		{
		protein->H[0][j]=protein->N[0][j]+bond[j]*CAT_Rbond_NH/b;
		}
	}
	for(j=0;j<3;j++) { bond[j]=protein->CA[0][j]-protein->N[0][j]; }
	b=norm_d(&(bond[0]),3);
	if(fabs(b-CAT_Rbond_CaN)>1e-12)
	{
		for(j=0;j<3;j++)
		{
		protein->CA[0][j]=protein->N[0][j]+bond[j]*CAT_Rbond_CaN/b;
		}
	}

	for(i=0;i<protein->n_res;i++)
	{
		CAT_add_peptide(protein,i,protein->phi[i],protein->psi[i]);
	}
}
void CAT_rescale( cat_prot *p)
{
	int i,j;
	double r;
	r=dist_d(p->N[0],p->CA[0],3);
	p->N[0][0]=p->CA[0][0] + CAT_Rbond_CaN/r*(p->N[0][0]-p->CA[0][0]);
	p->N[0][1]=p->CA[0][1] + CAT_Rbond_CaN/r*(p->N[0][1]-p->CA[0][1]);
	p->N[0][2]=p->CA[0][2] + CAT_Rbond_CaN/r*(p->N[0][2]-p->CA[0][2]);
	r=dist_d(p->H[0],p->N[0],3);
	p->H[0][0]=p->N[0][0] + CAT_Rbond_NH/r*(p->H[0][0]-p->N[0][0]);
	p->H[0][1]=p->N[0][1] + CAT_Rbond_NH/r*(p->H[0][1]-p->N[0][1]);
	p->H[0][2]=p->N[0][2] + CAT_Rbond_NH/r*(p->H[0][2]-p->N[0][2]);

	//MAH
	p->phi[0]=calc_dihedralf_angle(p->H[0],p->N[0],p->CA[0],p->C[0])-M_PI;
	p->psi[0]=calc_dihedralf_angle(p->N[0],p->CA[0],p->C[0],p->N[1]);
	for(i=1;i<p->n_res-1;i++)
	{
		p->phi[i]=calc_dihedralf_angle(p->C[i-1],p->N[i],p->CA[i],p->C[i]);
		p->psi[i]=calc_dihedralf_angle(p->N[i],p->CA[i],p->C[i],p->N[i+1]);
	}
	p->phi[i]=calc_dihedralf_angle(p->C[i-1],p->N[i],p->CA[i],p->C[i]);
	//FINE AGGIORNAMENTO DIEDRI
	double angle_NCaC[p->n_res];
	for(i=0;i<p->n_res-1;i++){
		angle_NCaC[i]=angle_ABC(p->N[i],p->CA[i],p->C[i]);
	}
	for(i=0;i<p->n_res;i++)
	{
		fprintf(stderr,"%7.5f %7.5f\t",p->phi[i],p->psi[i]);
	}
	for(i=1;i<p->n_res;i++){
		CAT_add_peptide(p,i,p->phi[i],p->psi[i],angle_NCaC[i]);
	}
	p->psi[0]=calc_dihedralf_angle(p->N[0],p->CA[0],p->C[0],p->N[1]);
	for(i=1;i<p->n_res-1;i++)
	{
		p->phi[i]=calc_dihedralf_angle(p->C[i-1],p->N[i],p->CA[i],p->C[i]);
		p->psi[i]=calc_dihedralf_angle(p->N[i],p->CA[i],p->C[i],p->N[i+1]);
	}
	p->phi[i]=calc_dihedralf_angle(p->C[i-1],p->N[i],p->CA[i],p->C[i]);
	for(i=0;i<p->n_res;i++)
	{
		fprintf(stderr,"%7.5f %7.5f\t",p->phi[i],p->psi[i]);
	}
}
*/

/**
 * @brief Function correct backbone distortions
 *
 * Tested OK. 03.08.16. (A shrinked Caterpillar protein is reinflated to itself)
 * 
 * @param[in,out]   *protein        cat_prot structure to be remodeled
 *
 * @return \c void
 */
void CAT_rescale( cat_prot *protein)
{
	int i,j,k,l;
	double bond[3];
	double b;
	double r_NH=CAT_Rbond_NH;
	double r_CO=CAT_Rbond_CO;
	double r_CN=CAT_Rbond_CN;
	double r_CCa=CAT_Rbond_CCa;
	double r_CaN=CAT_Rbond_CaN;

	for(i=0;i<protein->n_res;i++)
	{
		//Hydrogens
		for(j=0;j<3;j++) { bond[j]=protein->H[i][j]-protein->N[i][j]; }
		b=norm_d(&(bond[0]),3);
		if(fabs(b-CAT_Rbond_NH)>1e-10)
		{
			for(j=0;j<3;j++) { protein->H[i][j]=protein->N[i][j]+bond[j]*CAT_Rbond_NH/b;}
		}
		//Oxygens
		for(j=0;j<3;j++){bond[j]=protein->O[i][j]-protein->C[i][j];}
		b=norm_d(&(bond[0]),3);
		//if(fabs(b-CAT_Rbond_CO)>1e-10)
		if(fabs(b-r_CO)>1e-10)
		{
			for(j=0;j<3;j++) { protein->O[i][j]=protein->C[i][j]+bond[j]*r_CO/b;}
		}
		//Backbone atoms. N, Ca,C.
		//Nytrogens.
		for(j=0;j<3;j++){bond[j]=protein->CA[i][j]-protein->N[i][j];}
		b=norm_d(&(bond[0]),3);
		if(fabs(b-CAT_Rbond_CaN)>1e-10)
		{
			for(j=0;j<3;j++) { bond[j]*=(1.0-CAT_Rbond_CaN/b); }
			for(j=0;j<3;j++)
			{
				protein->CA[i][j]-=bond[j];
				protein->	C[i][j]-=bond[j];
				protein->	O[i][j]-=bond[j];
				protein->	H[i][j]-=bond[j];
				if (protein->CB!=NULL){
					protein->CB[i][j]-=bond[j];
				}
			}
			for(k=i+1;k<protein->n_res;k++)
			{
				for(j=0;j<3;j++)
				{
					protein->	N[k][j]-=bond[j];
					protein->CA[k][j]-=bond[j];
					protein->	C[k][j]-=bond[j];
					protein->	O[k][j]-=bond[j];
					protein->	H[k][j]-=bond[j];
					if(protein->CB!=NULL) {
						protein->CB[k][j]-=bond[j];
					}
				}
			}
		}
		if(i>0)
		{
			for(j=0;j<3;j++){bond[j]=protein->C[i-1][j]-protein->N[i][j];}
			b=norm_d(&(bond[0]),3);
			if(fabs(b-CAT_Rbond_CN)>1e-10)
			{
				for(j=0;j<3;j++) { bond[j]*=(1.0-CAT_Rbond_CN/b); }
				for(k=i;k<protein->n_res;k++)
				{
					for(j=0;j<3;j++)
					{
						protein->	N[k][j]	+=bond[j];
						protein->CA[k][j]	+=bond[j];
						protein->	C[k][j]	+=bond[j];
						protein->	O[k][j]	+=bond[j];
						protein->	H[k][j]	+=bond[j];
						if(protein->CB!=NULL){
							protein->CB[k][j]	+=bond[j];
						}
					}
				}
			}
		}
		//Carbon-alpha
		for(j=0;j<3;j++){bond[j]=protein->C[i][j]-protein->CA[i][j];}
		b=norm_d(&(bond[0]),3);
		if(fabs(b-CAT_Rbond_CCa)>1e-10)
		{
			for(j=0;j<3;j++) { bond[j]*=(1.0-CAT_Rbond_CCa/b); }
			for(j=0;j<3;j++)
			{
				protein->	C[i][j]-=bond[j];
				protein->	O[i][j]-=bond[j];
				protein->	H[i][j]-=bond[j];
				if(protein->CB!=NULL){
					protein->CB[i][j]-=bond[j];
				}
			}

			for(k=i+1;k<protein->n_res;k++)
			{
				for(j=0;j<3;j++)
				{
					protein->	N[k][j]	-=bond[j];
					protein->CA[k][j]	-=bond[j];
					protein->	C[k][j]	-=bond[j];
					protein->	O[k][j]	-=bond[j];
					protein->	H[k][j]	-=bond[j];
					if(protein->CB!=NULL) {
						protein->CB[k][j]	-=bond[j];
					}
				}
			}
		}
	}
	//Backbone angles.
	//Tested OK 02.11.16 (prova_rescale_angles.c)
	/*
	double dalpha;
	double alpha;
	double theta;
	double gamma;
	double bond_0[3];
	double bond_2[3];
	double axis[3];
	double ref_axis[3];
	double aNCaC = M_PI - CAT_angle_NCaC;
	double aCaCN = M_PI - CAT_angle_CaCN;
	double aCNCa = M_PI - CAT_angle_CNCa;
	for(i=0;i<protein->n_res;i++)
	{
		//--- CA
		//for(j=0;j<3;j++)
		//{
		//	bond[j]		=protein->C[i][j]-protein->CA[i][j];
		//	bond_0[j] =protein->CA[i][j]-protein->N[i][j];
		//}
		//alpha=acos(scal_d(bond,bond_0,3)/(r_CCa*r_CaN));
		//dalpha=alpha-aNCaC;
		//vecprod_d(bond,bond_0,axis);
		//k=(i+1)*protein->n_atom_per_res;
		////I have to do this to avoid rotating the hydrogen..
		//if(fabs(dalpha)>1e-12)
		//{
		//	if(i<protein->n_res-1)
		//	{
		//		rotate(&protein->coord[k],protein->n_atoms-k,protein->CA[i],axis,dalpha);
		//	}
		//	rotate(&protein->C[i],1,protein->CA[i],axis,dalpha);
		//	rotate(&protein->O[i],1,protein->CA[i],axis,dalpha);
		//}
		//--- C
		if(i<protein->n_res-1)
		{
			for(j=0;j<3;j++)
			{
				bond_0[j] =protein->C[i][j]-protein->CA[i][j];
				bond[j]		=protein->N[i+1][j]-protein->C[i][j];
				bond_2[j] =protein->O[i][j]-protein->C[i][j];
			}
			alpha=acos(scal_d(bond,bond_0,3)/(r_CCa*r_CN));
			gamma=acos(scal_d(bond_2,bond,3)/(r_CO*r_CN));
			theta=acos(scal_d(bond_0,bond_2,3)/(r_CO*r_CCa));
			if(gamma<theta) { dalpha=(alpha+aCaCN);}
			else  { dalpha=aCaCN-alpha;}
			if(fabs(dalpha)>1e-12)
			{
				vecprod_d(bond_2,bond_0,axis);
				k=(i+1)*protein->n_atom_per_res;
				rotate(&protein->coord[k],protein->n_atoms-k,protein->C[i],axis,dalpha);
			}
		}
		//--- N
		if(i>0)
		{
			for(j=0;j<3;j++)
			{
				bond_0[j] =protein->N[i][j] -protein->C[i-1][j];
				bond[j]		=protein->CA[i][j]-protein->N[i][j];
				bond_2[j] =protein->H[i][j]-protein->N[i][j];
			}
			alpha=acos(scal_d(bond,bond_0,3)/(r_CaN*r_CN));
			gamma=acos(scal_d(bond_2,bond,3)/(r_NH*r_CaN));
			theta=acos(scal_d(bond_0,bond_2,3)/(r_CN*r_NH));
			if(gamma<theta) { dalpha=(alpha+aCNCa);}
			else  { dalpha=aCNCa-alpha;}
			if(fabs(dalpha)>1e-12)
			{
				vecprod_d(bond_2,bond_0,axis);
				k=(i+1)*protein->n_atom_per_res;
				//I have to do this to avoid rotating the hydrogen..
				if(i<protein->n_res-1)
				{
					rotate(&protein->coord[k],protein->n_atoms-k,protein->N[i],axis,dalpha);
				}
				rotate(&protein->CA[i],1,protein->N[i],axis,dalpha);
				rotate(&protein->C[i],	1,protein->N[i],axis,dalpha);
				rotate(&protein->O[i],	1,protein->N[i],axis,dalpha);
			}
		}
	}
*/
 	CAT_prot_dihedrals (protein);
}

/*void CAT_add_peptide ( cat_prot * protein, int I, double phi, double psi, double angle_NCaC)
{
	int i,j;
	double a,b,c,cn;
	double alpha,beta,gamma,theta,theta_n;
	double norm;
	double z_old[3];
	double x[3],y[3],z[3]; //local basis
	double **X_loc=d2t(7,3);
	double **X_lab=d2t(7,3);
	double ax2[3];
	double P[3][3];
	double Rot[3][3];
	double M[3][3];
	double rq[4];

	protein->phi[I]=phi;
	protein->psi[I]=psi;
	//Atom position in the DH local base. origin on CA. reference dihedrals = pi.
	//C
	theta=CAT_angle_NCaC-M_PI_2;
	//PEZZOTTO
	//angle_NCaC=angle_ABC(protein->N[I],protein->CA[I],protein->C[I]);
	//theta=angle_NCaC-M_PI_2;
	//
	X_loc[0][0]=0;
	X_loc[0][1]=CAT_Rbond_CCa*cos(theta);
	X_loc[0][2]=CAT_Rbond_CCa*sin(theta);
	ax2[0]=0;
	ax2[1]=cos(theta);
	ax2[2]=sin(theta);
	//O
	a=CAT_Rbond_CCa;
	b=CAT_Rbond_CO;
	c=sqrt(a*a+b*b-2*a*b*cos(CAT_angle_CaCO));
	beta=asin(b/c*sin(CAT_angle_CaCO));
	theta-=beta;
	X_loc[1][0]=0;
	X_loc[1][1]= c*cos(theta);
	X_loc[1][2]= c*sin(theta);
	//N
	a=CAT_Rbond_CCa;
	b=CAT_Rbond_CN;
	c=sqrt(a*a+b*b-2*a*b*cos(CAT_angle_CaCN));
	//c=sqrt(a*a+b*b-2*a*b*cos(angle_CaCN));
	cn=c;
	beta=asin(b/c*sin(CAT_angle_CaCN));
	theta= CAT_angle_NCaC-M_PI_2+beta; //mind the sign
	//beta=asin(b/c*sin(CAT_angle_CaCN));
	//theta= angle_NCaC-M_PI_2+beta; //mind the sign
	theta_n=theta;
	X_loc[2][0]=0;
	X_loc[2][1]= c*cos(theta);
	X_loc[2][2]= c*sin(theta);
	//H
	alpha=asin(a/cn*sin(CAT_angle_CaCN));
	//alpha=asin(a/cn*sin(angle_CaCN));
	gamma=CAT_angle_CNH-alpha;
	//change triangle to CaNH
	a=cn;
	b=CAT_Rbond_NH;
	c=sqrt(a*a+b*b-2*a*b*cos(gamma));
	beta=asin(b/c*sin(gamma));
	theta=theta_n+beta;
	X_loc[3][0]=0;
	X_loc[3][1]= c*cos(theta);
	X_loc[3][2]= c*sin(theta);
	//CA
	alpha=asin(CAT_Rbond_CCa/cn*sin(CAT_angle_CaCN));
	gamma=alpha+CAT_angle_CNCa;
	a=cn;
	b=CAT_Rbond_CaN;
	c=sqrt(a*a+b*b-2*a*b*cos(gamma));
	beta=asin(b/c*sin(gamma));
	theta=theta_n-beta;
	X_loc[4][0]=0;
	X_loc[4][1]= c*cos(theta);
	X_loc[4][2]= c*sin(theta);
	//a rigor di logica..per simmetria...
	alpha=M_PI-gamma-beta;
	gamma=CAT_angle_NCaC+alpha;
	//PEZZOTTO
	//if(I<protein->n_res-1){
	//	angle_NCaC=angle_ABC(protein->N[I+1],protein->CA[I+1],protein->C[I+1]);
	//}else{
//		angle_NCaC=CAT_angle_NCaC;
	//}
	//--
	//--
	a=c;
	b=CAT_Rbond_CCa;
	c=sqrt(a*a+b*b-2*a*b*cos(gamma));
	beta=asin(b/c*sin(gamma));
	theta+=beta;
	X_loc[5][0]=0;
	X_loc[5][1]=c*cos(theta);
	X_loc[5][2]=c*sin(theta);
	//O
	alpha=M_PI-gamma-beta;
	gamma=CAT_angle_CaCO-alpha;
	a=c;
	b=CAT_Rbond_CO;
	c=sqrt(a*a+b*b-2*a*b*cos(gamma));
	beta=asin(b/c*sin(gamma));
	theta+=beta;
	X_loc[6][0]=0;
	X_loc[6][1]=c*cos(theta);
	X_loc[6][2]=c*sin(theta);

	double zero[3]={0,0,0};
	double axis[3]={-1,0,0};
	rotate(X_loc,7,zero,axis,CAT_angle_NCaC-angle_NCaC);
	//Rotation matrix moving X_loc to have angles phi,psi given in input
	//quaternion from multiplication of the two rotations.
	phi=0.5*phi;
	psi=0.5*psi;
	rq[0]= sin(phi)*sin(psi) - ax2[2]*cos(psi)*cos(phi);
	rq[1]= -cos(psi)*( sin(phi)*ax2[0] + cos(phi)*ax2[1]);
	rq[2]= -cos(psi)*( sin(phi)*ax2[1] - cos(phi)*ax2[0]);
	rq[3]= -cos(phi)*sin(psi)-cos(psi)*sin(phi)*ax2[2];
	//rotation matrix
	//first col
	Rot[0][0]=1-2*(rq[2]*rq[2]+rq[3]*rq[3]);
	Rot[0][1]=2*(rq[2]*rq[1]+rq[0]*rq[3]);
	Rot[0][2]=2*(rq[3]*rq[1]-rq[0]*rq[2]);
	//second col
	Rot[1][0]=2*(rq[1]*rq[2]-rq[0]*rq[3]);
	Rot[1][1]=1-2*(rq[3]*rq[3]+rq[1]*rq[1]);
	Rot[1][2]=2*(rq[3]*rq[2]+rq[0]*rq[1]);
	//third col
	Rot[2][0]=2*(rq[1]*rq[3]+rq[0]*rq[2]);
	Rot[2][1]=2*(rq[3]*rq[2]-rq[0]*rq[1]);
	Rot[2][2]=1-2*(rq[2]*rq[2]+rq[1]*rq[1]);
	//Matrix to change ref. frame
	if(I>0)
	{
		z_old[0]=(protein->N[I][0]-protein->C[I-1][0]);
		z_old[1]=(protein->N[I][1]-protein->C[I-1][1]);
		z_old[2]=(protein->N[I][2]-protein->C[I-1][2]);
	}
	else
	{
		z_old[0]=0;
		z_old[1]=0;
		z_old[2]=1;
	}
	z[0]=(protein->CA[I][0]-protein->N[I][0]);
	z[1]=(protein->CA[I][1]-protein->N[I][1]);
	z[2]=(protein->CA[I][2]-protein->N[I][2]);
	norm=sqrt(z[0]*z[0]+z[1]*z[1]+z[2]*z[2]);
	z[0]/=norm; z[1]/=norm; z[2]/=norm;
	//
  x[0] = z_old[1] * z[2] - z_old[2] * z[1];
  x[1] = z_old[2] * z[0] - z_old[0] * z[2];
  x[2] = z_old[0] * z[1] - z_old[1] * z[0];
	norm=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
	x[0]/=norm; x[1]/=norm; x[2]/=norm;
	//
  y[0] = z[1] * x[2] - z[2] * x[1];
  y[1] = z[2] * x[0] - z[0] * x[2];
  y[2] = z[0] * x[1] - z[1] * x[0];
	//[col][row] again
	P[0][0]=x[0]; P[0][1]=x[1]; P[0][2]=x[2];
	P[1][0]=y[0]; P[1][1]=y[1]; P[1][2]=y[2];
	P[2][0]=z[0]; P[2][1]=z[1]; P[2][2]=z[2];
	//Putting all together
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			//having transposed the matrix order is a bit confusing
			M[j][i] =P[0][i]*Rot[j][0];
			M[j][i]+=P[1][i]*Rot[j][1];
			M[j][i]+=P[2][i]*Rot[j][2];
		}
	}
	for(j=0;j<7;j++)
	{
		X_lab[j][0]=0;
		X_lab[j][1]=0;
		X_lab[j][2]=0;
		for(i=0;i<3;i++)
		{
			X_lab[j][0]+=M[i][0]*X_loc[j][i];
			X_lab[j][1]+=M[i][1]*X_loc[j][i];
			X_lab[j][2]+=M[i][2]*X_loc[j][i];
		}
		//shift
		X_lab[j][0]+=protein->CA[I][0];
		X_lab[j][1]+=protein->CA[I][1];
		X_lab[j][2]+=protein->CA[I][2];
	}
	//rotate to maintain bending
	double v[3],w[3];
	double n[3];
	for(j=0;j<3;j++){
		v[j]=protein->CA[I][j]-protein->N[I][j];
		w[j]=protein->C[I][j] -protein->CA[I][j];
	}
	vecprod_d(v,w,n);
	//rotate(X_lab,7,protein->CA[I],n,CAT_angle_NCaC-angle_NCaC);
	//qui ci vuole una condizione per vedere se sto andando fuori dai boundaries
	for(j=0;j<3;j++)
	{
		protein->C   [I][j]= X_lab[0][j];
		protein->O   [I][j]= X_lab[1][j];
		if(I<protein->n_res-1)
		{
			protein->N [I+1][j]= X_lab[2][j];
			protein->H [I+1][j]= X_lab[3][j];
			protein->CA[I+1][j]= X_lab[4][j];
			protein->C [I+1][j]= X_lab[5][j];
			protein->O[I+1][j] = X_lab[6][j];
		}
	}
	if(protein->CB!=NULL)
	{
		CAT_insert_cbeta(protein,I,CAT_Cb_AZIMUTH,CAT_Rbond_CCb);
	}
	free_d2t(X_lab);
	free_d2t(X_loc);
}
*/

/**
 * @brief Function modifie residue configuration base on dihedrals
 *
 * Function reconstruct given residue based on dihedral angles in cat_prot structure.
 * Function does not change dihedrals in protein phi/psi have to be recalculated.
 * Function do not check if change cause breaking of protein backbone conectivity!
 * Tested OK. LT 02.08.16
 * 
 * @param[in,out]   *p              cat_prot structure to be remodeled
 * @param[in]        I              Residue which atom positions are recalculated
 * @param[in]        phi            PHI angle (in radians) of given residue
 * @param[in]        alpha          Deviation of 
 * @param[in]        psi            PSI angle (in radians) of given residue
 *
 * @return err code
 */
int CAT_add_peptide ( cat_prot *p, int I, double phi, double alpha, double psi )
{
	int error;
	double q1[4]={0,0,0,0};
	double q2[4]={0,0,0,0};
	double q3[4]={0,0,0,0};
	double qmul[4]={0,0,0,0};
	double pep[3*6];
	//matrices
	double Pm[9];
	double rot[9];
	double Rt[16];
	gsl_matrix_view 			pep_v  	=gsl_matrix_view_array(pep,3,6);
	gsl_matrix_view 			Pm_v		=gsl_matrix_view_array(Pm,3,3);
	gsl_matrix_view 			rot_v 	=gsl_matrix_view_array(rot,3,3);
	gsl_matrix_view				Rt_v 		=gsl_matrix_view_array(Rt,4,4);
    gsl_matrix_set_all(&Rt_v.matrix,0.0); // Explicit initialization of Rt(rotation matrix) ... if we want to make valgrind happy
	gsl_matrix_view				M_v  		=gsl_matrix_submatrix (&Rt_v.matrix,0,0,3,3);
	//quaternions
	gsl_vector_view q1_v 	=gsl_vector_view_array(q1,4);
	gsl_vector_view q2_v 	=gsl_vector_view_array(q2,4);
	gsl_vector_view q3_v	=gsl_vector_view_array(q3,4);
	gsl_vector_view qmul_v	=gsl_vector_view_array(qmul,4);
	double axPhi[3]={0,0,1};
	double axAlpha[3]={-1,0,0};
	double axPsi[3];
	phi-=M_PI; psi-=M_PI;
	gsl_vector_view axAlpha_v =gsl_vector_view_array(axAlpha,3);
	gsl_vector_view axPsi_v		=gsl_vector_view_array(axPsi,3);
	gsl_vector_view axPhi_v		=gsl_vector_view_array(axPhi,3);
	//get ideal peptide
	build_peptide (&pep_v.matrix); //returns CA-C-O-N-H-CA(I+1)
	//Phi
	quaternion_build(&q1_v.vector,&axPhi_v.vector,phi);
	//Alpha
	quaternion_build(&q2_v.vector,&axAlpha_v.vector,alpha);
	//Psi
	gsl_vector_view v_v =gsl_matrix_column(&pep_v.matrix,1);
	gsl_vector_memcpy(&axPsi_v.vector,&v_v.vector);
	quaternion_build(&q3_v.vector,&axPsi_v.vector,psi);
	//multiply quaternions
	quaternion_mult (&qmul_v.vector,&q2_v.vector,&q3_v.vector);
	gsl_vector_memcpy(&q2_v.vector,&qmul_v.vector);
	quaternion_mult (&qmul_v.vector,&q1_v.vector,&q2_v.vector);
	//rotation matrix
	rotation3D_build(&rot_v.matrix,&qmul_v.vector);
	//Tranformation matrix to lab ref. system
	double z_old[3];
	double z_lab[3];
	double x_lab[3];
	double y_lab[3];
	if(I>0){
		for(int i=0;i<3;i++) { z_old[i]=p->N[I][i]-p->C[I-1][i];}
	}
	else //QUESTA VA RISOLTA!!
	{
		//reconstruct z_old from the projection of the ideal bond on 
		//the bonds NH and NCa;
		double b1[3],b2[3];
		for(int i=0;i<3;i++){
			b1[i]=p->H[0][i]-p->N[0][i];
			b2[i]=p->CA[0][i]-p->N[0][i];
		}
		normalize_d(b1,3);
		normalize_d(b2,3);
		for(int i=0;i<3;i++){
			z_old[i]=-cos(CAT_angle_CNH)*b1[i];
			z_old[i]-=cos(CAT_angle_CNCa)*b2[i];
		}
	}
	for(int i=0;i<3;i++) { z_lab[i]=p->CA[I][i]-p->N[I][i];}
	normalize_d(z_lab,3);
	vecprod_d(z_old,z_lab,x_lab);
	normalize_d(x_lab,3);
	normalize_d(z_lab,3);
	vecprod_d(z_lab,x_lab,y_lab);
	normalize_d(y_lab,3);
	for(int i=0;i<3;i++) {
		gsl_matrix_set(&Pm_v.matrix,i,0,x_lab[i]);
		gsl_matrix_set(&Pm_v.matrix,i,1,y_lab[i]);
		gsl_matrix_set(&Pm_v.matrix,i,2,z_lab[i]);
	}
	//compose with rotation
	error=gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,&Pm_v.matrix,&rot_v.matrix,0,&M_v.matrix);
	//add translation
	gsl_matrix_set(&Rt_v.matrix,0,3,p->CA[I][0]);
	gsl_matrix_set(&Rt_v.matrix,1,3,p->CA[I][1]);
	gsl_matrix_set(&Rt_v.matrix,2,3,p->CA[I][2]);
	gsl_matrix_set(&Rt_v.matrix,3,3,1);
	//perform rototranslation
	double gen_coord		[4*6];
	double gen_coord_rot[4*6];
	gsl_matrix_view gc_v 	= gsl_matrix_view_array(gen_coord,4,6);
	gsl_matrix_view gcr_v = gsl_matrix_view_array(gen_coord_rot,4,6);
	gsl_matrix_set_all(&gc_v.matrix,1.0);
    gsl_matrix_set_all(&gcr_v.matrix,0.0); // Explicit initialization of gcr_v(result matrix) ... if we want to make valgrind happy
	gsl_matrix_view c_v= gsl_matrix_submatrix(&gc_v.matrix,0,0,3,6);
	gsl_matrix_memcpy(&c_v.matrix,&pep_v.matrix);
	error=gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,&Rt_v.matrix,&gc_v.matrix,0.,&gcr_v.matrix);
	//get coordinates
	for(int i=0;i<3;i++){
		p->C[I][i]		=gsl_matrix_get(&gcr_v.matrix,i,1);
		p->O[I][i]		=gsl_matrix_get(&gcr_v.matrix,i,2);
		if(I+1<p->n_res){
			p->N[I+1][i]	=gsl_matrix_get(&gcr_v.matrix,i,3);
			p->H[I+1][i]	=gsl_matrix_get(&gcr_v.matrix,i,4);
			p->CA[I+1][i]	=gsl_matrix_get(&gcr_v.matrix,i,5);
		}
	}
	return error;
}

/**
 * @brief Function recalculate protein dihedral angles
 *
 * Function iterate through whole protein and recalculate p->phi and p->psi.
 * 
 * @param[in,out]   *p        Protein which dihedrals are recalculated
 *
 * @return \c void
 */
void CAT_prot_dihedrals (cat_prot *protein)
{
	int i;
	for(i=0;i<protein->n_res;i++)
	{
		if(i>0)
		{
			protein->phi[i]=calc_dihedralf_angle(protein->C[i-1],protein->N[i],protein->CA[i],protein->C[i]);
		}
		if(i<protein->n_res-1)
		{
			protein->psi[i]=calc_dihedralf_angle(protein->N[i],protein->CA[i],protein->C[i],protein->N[i+1]);
		}
	}
}

/**
 * @brief Function calculate phi angle of residue in protein
 *
 *
 * @param[in]      *p        Protein from where angle is calculated
 * @param[in]       c        Index of atom which dihedral is calculated
 *
 * @return phi dihedral angle
 */
double compute_phi(cat_prot *p,int c)
{
	int k;
	double cs, sn;
	double v_1[3],v_2[3],w_1[3],w_2[3],w_3[3];
	double axis[3];
	//phi
	if(c>1)
	{
		for(k=0;k<3;k++)
		{
			v_1 [k] = p->C [c-1][k] -	p->N [c][k];
		}
	}
	else
	{
		for(k=0;k<3;k++)
		{
			v_1 [k] = p->H [c][k] -	p->N [c][k];
		}
	}
	for(k=0;k<3;k++)
	{
		v_2 [k] = p->C [c][k]		-	p->CA[c][k];
		axis[k] = p->CA[c][k] 	- p->N [c][k];
	}
	vecprod_d(v_1,axis,w_1);normalize_d(w_1,3);
	vecprod_d(v_2,axis,w_2);normalize_d(w_2,3);
	cs=scal_d(w_1,w_2,3);
	vecprod_d(w_1,w_2,w_3);
	sn=norm_d(w_3,3)*GSL_SIGN(scal_d(w_3,axis,3));
	//p->phi[c]=atan2(sn,cs);
	return atan2(sn,cs);
}

/**
 * @brief Function calculate psi angle of residue in protein
 *
 *
 * @param[in]      *p        Protein from where angle is calculated
 * @param[in]       c        Index of atom which dihedral is calculated
 *
 * @return psi dihedral angle
 */
double compute_psi(cat_prot *p,int c)
{
	int k;
	double cs, sn;
	double v_1[3],v_2[3],w_1[3],w_2[3],w_3[3];
	double axis[3];
	//psi
	if(c<p->n_res-1)
	{
		for(k=0;k<3;k++)
		{
			v_2 [k] = p->N [c+1][k]	-	p->C [c][k];
		}
	}
	else
	{
		for(k=0;k<3;k++)
		{
			v_2 [k] = p->O [c][k]	-	p->C [c][k];
		}
	}
	for(k=0;k<3;k++)
	{
		v_1 [k] = p->N [c][k] 	-	p->CA[c][k];
		axis[k] = p->C [c][k] 	- p->CA[c][k];
	}
	vecprod_d(v_1,axis,w_1);normalize_d(w_1,3);
	vecprod_d(v_2,axis,w_2);normalize_d(w_2,3);
	cs=scal_d(w_1,w_2,3);
	vecprod_d(w_1,w_2,w_3);
	sn=norm_d(w_3,3)*GSL_SIGN(scal_d(w_3,axis,3));
	//p->psi[c]=atan2(sn,cs);
	return atan2(sn,cs);
}

/**
 * @brief Function recalculate protein dihedral angles
 *
 * Function iterate through whole protein and recalculate p->phi and p->psi.
 * 
 * @param[in,out]   *p        Protein which dihedrals are recalculated
 *
 * @return \c void
 */
void compute_dihedrals(cat_prot *p)
{
	int c,k;
	double cs, sn;
	double v_1[3],v_2[3],w_1[3],w_2[3],w_3[3];
	double axis[3];
	//phi
	for(c=1;c<p->n_res;c++)
	{
		for(k=0;k<3;k++) 
		{ 
			v_1 [k] = p->C [c-1][k] -	p->N [c][k]; 
			v_2 [k] = p->C [c][k]		-	p->CA[c][k]; 
			axis[k] = p->CA[c][k] 	- p->N [c][k];
		}
		vecprod_d(v_1,axis,w_1);normalize_d(w_1,3);
		vecprod_d(v_2,axis,w_2);normalize_d(w_2,3);
		cs=scal_d(w_1,w_2,3);
		vecprod_d(w_1,w_2,w_3); 
		sn=norm_d(w_3,3)*GSL_SIGN(scal_d(w_3,axis,3));
		p->phi[c]=atan2(sn,cs);
	}
	//psi
	for(c=0;c<p->n_res-1;c++)
	{
		for(k=0;k<3;k++) 
		{ 
			v_1 [k] = p->N [c][k] 	-	p->CA[c][k]; 
			v_2 [k] = p->N [c+1][k]	-	p->C [c][k]; 
			axis[k] = p->C [c][k] 	- p->CA[c][k];
		}
		vecprod_d(v_1,axis,w_1);normalize_d(w_1,3);
		vecprod_d(v_2,axis,w_2);normalize_d(w_2,3);
		cs=scal_d(w_1,w_2,3);
		vecprod_d(w_1,w_2,w_3); 
		sn=norm_d(w_3,3)*GSL_SIGN(scal_d(w_3,axis,3));
		p->psi[c]=atan2(sn,cs);
	}
}

/**
 * @brief Function calculate dihedral from 4 atoms
 *
 * This function was taken from LAMMPS.
 *
 * @param[in]      *atom_1        Coordinates of atom1 in 3D
 * @param[in]      *atom_2        Coordinates of atom2 in 3D
 * @param[in]      *atom_3        Coordinates of atom3 in 3D
 * @param[in]      *atom_4        Coordinates of atom4 in 3D
 *
 * @return dihedral angle between atoms
 */
double calc_dihedralf_angle(double *atom_1, double *atom_2, double *atom_3, double *atom_4)
{
	double vb1x=0,vb1y=0,vb1z=0,vb2xm=0,vb2ym=0,vb2zm=0,vb3x=0,vb3y=0,vb3z=0;
	double ax=0,ay=0,az=0,bx=0,by=0,bz=0,rasq=0,rbsq=0,rab;
	double c=0,s=0;
	double rgsq,rg;
	double dihedral=0;
	char err_msg[1024];
	//double sinphi=0,rgsq=0;
	// 1st bond
	vb1x=atom_1[0]-atom_2[0];
	vb1y=atom_1[1]-atom_2[1];
	vb1z=atom_1[2]-atom_2[2];
	// 2nd bond (opposite direction -- axis)
	vb2xm=atom_2[0]-atom_3[0];
	vb2ym=atom_2[1]-atom_3[1];
	vb2zm=atom_2[2]-atom_3[2];
	// 3rd bond
	vb3x=atom_4[0]-atom_3[0];
	vb3y=atom_4[1]-atom_3[1];
	vb3z=atom_4[2]-atom_3[2];
	// c,s calculation
	ax = vb1y*vb2zm - vb1z*vb2ym;
	ay = vb1z*vb2xm - vb1x*vb2zm;
	az = vb1x*vb2ym - vb1y*vb2xm;

	bx = vb3y*vb2zm - vb3z*vb2ym;
	by = vb3z*vb2xm - vb3x*vb2zm;
	bz = vb3x*vb2ym - vb3y*vb2xm;

	rasq = ax*ax + ay*ay + az*az;
	rbsq = bx*bx + by*by + bz*bz;
	rgsq = vb2xm*vb2xm + vb2ym*vb2ym + vb2zm*vb2zm;
	rg = sqrt(rgsq);

	//ra2inv = 1.0/rasq;
	//rb2inv = 1.0/rbsq;
	//rabinv = sqrt(ra2inv*rb2inv);
	rab=sqrt(rasq*rbsq);

	c= (ax*bx + ay*by + az*bz)/rab;
	s= rg*(ax*vb3x + ay*vb3y + az*vb3z)/rab;


	//if(c < -1) c=-0.999999999999;
	//if(c > 1) c=0.9999999999999; // The DBL_EPSILON is there to make sure that the cosphi is always a bit smaller than one otherwise acos returns a NaN
	dihedral=atan2(s,c);
	if(isnan(dihedral))
	{
		sprintf(err_msg,"%s:%d  cosphi=%g sinphi=%g\n",__FILE__,__LINE__,c,s);
        printf("DEBUG: d0 %g %g %g | %g %g %g | %g %g %g | %g %g %g\n", atom_1[0], atom_1[1], atom_1[2], atom_2[0], atom_2[1], atom_2[2], atom_3[0], atom_3[1], atom_3[2],atom_4[0], atom_4[1], atom_4[2]);
		failed(err_msg);
	}
	// Calculate dihedral angle
	//if( scalar(aXb, c) < 0.0 ) *phi = (2.0*M_PI) - *phi;
	return (dihedral);
}

/**
 * @brief Don not know ... Luca?
 * * 
 * @param[in,out]   *pep            ???????????????
 *
 * @return \c void
 */
void build_peptide ( gsl_matrix *pep)
{
	double a,b,c,cn;
	double alpha,beta,gamma,theta,theta_n;
	//Atom position in the DH local base. origin on CA. reference dihedrals = pi.
	//Ca
	gsl_matrix_set(pep,0,0,0);
	gsl_matrix_set(pep,1,0,0);
	gsl_matrix_set(pep,2,0,0);
	//C
	c=CAT_Rbond_CCa;
	theta=M_PI_2;
	gsl_matrix_set(pep,0,1,0);
	gsl_matrix_set(pep,1,1,c*cos(theta));
	gsl_matrix_set(pep,2,1,c*sin(theta));
	//O
	a=CAT_Rbond_CCa;
	b=CAT_Rbond_CO;
	c=sqrt(a*a+b*b-2*a*b*cos(CAT_angle_CaCO));
	beta=asin(b/c*sin(CAT_angle_CaCO));
	theta-=beta;
	gsl_matrix_set(pep,0,2,0);
	gsl_matrix_set(pep,1,2,c*cos(theta));
	gsl_matrix_set(pep,2,2,c*sin(theta));
	//N
	a=CAT_Rbond_CCa;
	b=CAT_Rbond_CN;
	c=sqrt(a*a+b*b-2*a*b*cos(CAT_angle_CaCN));
	cn=c;
	beta=asin(b/c*sin(CAT_angle_CaCN));
	theta= M_PI_2+beta; //mind the sign
	theta_n=theta;
	gsl_matrix_set(pep,0,3,0);
	gsl_matrix_set(pep,1,3,c*cos(theta));
	gsl_matrix_set(pep,2,3,c*sin(theta));
	//H
	alpha=asin(a/cn*sin(CAT_angle_CaCN));
	gamma=CAT_angle_CNH-alpha;
	//change triangle to CaNH
	a=cn;
	b=CAT_Rbond_NH;
	c=sqrt(a*a+b*b-2*a*b*cos(gamma));
	beta=asin(b/c*sin(gamma));
	theta=theta_n+beta;
	gsl_matrix_set(pep,0,4,0);
	gsl_matrix_set(pep,1,4,c*cos(theta));
	gsl_matrix_set(pep,2,4,c*sin(theta));
	//CA
	alpha=asin(CAT_Rbond_CCa/cn*sin(CAT_angle_CaCN));
	gamma=alpha+CAT_angle_CNCa;
	a=cn;
	b=CAT_Rbond_CaN;
	c=sqrt(a*a+b*b-2*a*b*cos(gamma));
	beta=asin(b/c*sin(gamma));
	theta=theta_n-beta;
	gsl_matrix_set(pep,0,5,0);
	gsl_matrix_set(pep,1,5,c*cos(theta));
	gsl_matrix_set(pep,2,5,c*sin(theta));
}

