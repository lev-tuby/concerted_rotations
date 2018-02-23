#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf.h>
#include <sys/types.h>
#include <unistd.h>


#include "./lib/Caterpillar_energies.h"
#include "./lib/my_memory.h"
#include "./lib/my_geom.h"
#include "./lib/geom_prop.h"
#include "./lib/Caterpillar_IO.h"
#include "./lib/CAT_moves.h"
#include "./lib/quaternions.h"
#include "./lib/histogram.h"

#define ACC 1
#define REJ 0

#define GLOBAL_SEED 100
#define NUMBER_OF_RESIDUES 20
#define SEQUENCE "AAAAAAAAAAAAAAAAAAAA"

// Values used to check backbone geometry ... maximal deviation in ANGLE BOND and OMEGA dihedral
#define ANGLE_MAX_DEVIATION 10e-9
#define BOND_MAX_DEVIATION  10e-9
#define OMEGA_MAX_DEVIATION 10e-9

#define CATENR_CaCa_SPRING 40

typedef struct {
    cat_prot *p;
    double E;
    double Ewat;
    double DE;
    double *Dcontacts;
}conf;

typedef struct {
	int step;
	gsl_rng *rng_r;
	double temp;
	double beta;
	conf old;
	conf new;
}mc_traj_data;


// Numerical stability control
int prot_rescale ( cat_prot *p);
void print_bond_errors(FILE *stream,cat_prot *p);
void print_joint_angles_errors(FILE *stream, cat_prot *p);
void print_omega_errors(FILE *stream, cat_prot *p);


mc_traj_data * mc_traj_data_alloc(int seed,double temp);
void mc_traj_data_free( mc_traj_data * mc_trj);
int acceptMove(mc_traj_data *mctrj);

void CAT_copy(cat_prot *dest, cat_prot *orig);
void Init_MC(mc_move_data **mvdt, mc_traj_data **mc_traj, char file_conf[1024],char file_pot[1024]);

    int global_warn=0;
    int global_max_warn=0;


int main(int argc, char *argv[])
{
    int
        maxIter=10e8;

    double
        sigma=0.1,
        *point;
    point = (double*)calloc(2, sizeof(double));

    char
        histogramName[1024];

    FILE
        *dihed_out=fopen("dihedrals.dat","w");

    mc_traj_data
        *mc_trj;

    histogram
        *psi_phi = histogram_init  (-M_PI, M_PI+0.000001, 180, 0.1);
    histogram_add_dimension(-M_PI, M_PI, 180+0.000001, psi_phi);

    mc_move_data
        *mc_mvdt;
    Init_MC(&mc_mvdt, &mc_trj, "unused",argv[1]);

    //MAIN LOOP
    for(int i=0;i<maxIter;i++)
    {
        // geometry check
        if(i>0 && i%1000==0)
        {
            fprintf(stdout,"iteration: %d\n",i);
            print_bond_errors(stdout,mc_trj->new.p);
            print_joint_angles_errors(stdout,mc_trj->new.p);
            print_omega_errors(stdout,mc_trj->new.p);
        }

        // rescale which fix numerical inaccuracies introduced in moves and recalculate system energy
        if(i%100000)
        {
            prot_rescale(mc_trj->old.p);
            prot_rescale(mc_trj->new.p);
        }

        // Histogram sampling
        if(i>1000000 && i%50==0)
        {
            for (size_t resi=1; resi < mc_trj->new.p->n_res-1; resi++)
            {
                point[0]=mc_trj->old.p->psi[resi];
                point[1]=mc_trj->old.p->phi[resi];
                histogram_add(point, psi_phi);
                fprintf(dihed_out,"%lf %lf  ",mc_trj->new.p->phi[resi],mc_trj->new.p->psi[resi]);
            }
            fprintf(dihed_out,"\n");
        }

        // Print out histogram and protein in pdb
        if(i%100000==0)
        {
            snprintf(histogramName, sizeof(histogramName), "sigma_%03.6lf_iteration_%08i.dat", sigma, i);
            histogram_print(histogramName, "matrix2d", psi_phi);
            CATIO_cat2pdb("prova.pdb","a","--",mc_trj->old.p,1); // save current configuratin
        }

        // MOVES
        if(i<5000) //first 5000 steps equil just with pivot move
        {
            CATMV_pivot(mc_mvdt,mc_trj->new.p,mc_trj->rng_r, mc_trj->new.p->n_res/2);
        } else if (gsl_rng_uniform(mc_trj->rng_r)>0.5) {
            CATMV_pivot(mc_mvdt,mc_trj->new.p,mc_trj->rng_r, mc_trj->new.p->n_res/5);
        } else {
            CATMV_concerted_rot(mc_mvdt,mc_trj->new.p,mc_trj->rng_r, sigma);
        }

        acceptMove(mc_trj); // copy configuration altered configuration in current configration
    } // END MAIN LOOP

    // Free used memory 
    mc_traj_data_free(mc_trj);
    CATMV_mc_move_data_free(mc_mvdt);
    histogram_free(psi_phi);
    return 0;
}



mc_traj_data * mc_traj_data_alloc(int seed,double temp)
{
    mc_traj_data
        *mctrj = (mc_traj_data*)malloc(sizeof(mc_traj_data));

    mctrj->step     = 0;
	mctrj->rng_r    = gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set(mctrj->rng_r, seed);
	mctrj->beta     = 1/temp; //Kb=1
	mctrj->temp     = temp;
	mctrj->old.p    = NULL;
	mctrj->old.E    = -1;
	mctrj->new.p    = NULL;
	mctrj->new.E    = -1;
	return mctrj;
}

void mc_traj_data_free( mc_traj_data * mctrj)
{
	if(mctrj!=NULL)
	{
		gsl_rng_free(mctrj->rng_r);
		free(mctrj);
	}
}


void CAT_copy(cat_prot *dest, cat_prot *orig)
{
	//Quite unefficient: I am copying the whole protein each time
	char err[128];
	if (NULL == dest )
	{
		failed("CAT_copy, NULL destination.\n");
	}
	if (NULL == orig )
	{
		failed("CAT_copy, NULL origin.\n");
	}
	if(dest->n_atoms != orig->n_atoms)
	{
		sprintf(err,
				"CAT_copy. Proteins have different numbers of atoms. dest: %lu. orig %lu",
				dest->n_atoms,orig->n_atoms);
		failed(err);
	}
	if(dest->n_atom_per_res != orig->n_atom_per_res)
	{
		sprintf(err,
				"CAT_copy. Proteins have different numbers of atoms per residue. dest: %lu. orig %lu",
				dest->n_atom_per_res,orig->n_atom_per_res);
		failed(err);
	}
	for(int i=0;i<orig->n_atoms;i++)
	{
		dest->coord[i][0]=orig->coord[i][0];
		dest->coord[i][1]=orig->coord[i][1];
		dest->coord[i][2]=orig->coord[i][2];
	}
	for(int i=0;i<orig->n_res;i++)
	{
		dest->phi[i]=orig->phi[i];
		dest->psi[i]=orig->psi[i];
		dest->residues[i]=orig->residues[i];
		dest->contacts[i]=orig->contacts[i];
	}
}


void Init_MC(mc_move_data **mvdt, mc_traj_data **mc_traj, char file_conf[1024],char file_pot[1024])
{
    int
        N_res=NUMBER_OF_RESIDUES,       // Number of residues in the protein
        N_atoms_per_res=5,              // Without CB N_atoms_per_res = 5 (N, C, O, H, CA), with CB N_atoms_per_res = 6 (N, C, O, H, CA, CB) !!!CB (untested)!!!
        seed = GLOBAL_SEED;             // Seed for initialization of random number generator

    double
        temp = 2.0,                     // Temperature in kT 
        C[3]={0.,0.,0.};                // Define position of first atom in peptide chain (used in initialization of extended peptide configuration)

    mc_move_data
        *mv = CATMV_mc_move_data_alloc(N_res,N_res);

    mc_traj_data
        *mct;
    mct         = mc_traj_data_alloc(seed, temp);
    mct->old.p  = CAT_prot_alloc(N_res, N_atoms_per_res);               // alocate structure for protein configuration
    mct->new.p  = CAT_prot_alloc(N_res,N_atoms_per_res);                // alocate structure for updates of configurations

    CAT_set_prot_linear(mct->old.p, C, 0.0);                            // create initial configuration of peptide backbone (extended in this case)
    CAT_set_residues_fasta(mct->old.p, N_res, SEQUENCE);                 // Initialize types of protein residues from fasta fromated string
    CAT_copy(mct->new.p,mct->old.p);                                    // now we copy configuration and residue types from old to new

	mct->old.E = 0.0;
	mct->old.Ewat=0;
	mct->old.Dcontacts=d1t(mct->old.p->n_res);


	mct->new.E=mct->old.E;
	mct->new.Ewat=mct->old.Ewat;
	mct->new.Dcontacts=d1t(mct->new.p->n_res);

	*mvdt=mv;
	*mc_traj=mct;
}

int acceptMove(mc_traj_data *mctrj)
{
    CAT_copy(mctrj->old.p,mctrj->new.p);
    mctrj->step++;
    return ACC;
}

int prot_rescale ( cat_prot *p)
{
    int
        error;

    double
        r,                  // used for length of bonds
        alpha[p->n_res],
        b1[3],
        b2[3];

	//ref=M_PI-CAT_angle_NCaC;

	for(int i=0;i<3;i++) {
		b1[i]=p->H[0][i]-p->N[0][i];
		b2[i]=p->N[0][i]-p->CA[0][i];
	}
	r=norm_d(b1,3);
	if(fabs(r-CAT_Rbond_NH)>1e-10) {
		for(int i=0;i<3;i++) {
			p->H[0][i]=p->N[0][i]+b1[i]*CAT_Rbond_NH/r;
		}
	}
	r=norm_d(b2,3);
	if(fabs(r-CAT_Rbond_CaN)>1e-10) {
		for(int i=0;i<3;i++){
			p->H[0][i]+=b2[i]*(CAT_Rbond_CaN/r-1.0);
			p->N[0][i]+=b2[i]*(CAT_Rbond_CaN/r-1.0);
		}
	}
	p->phi[0]=calc_dihedralf_angle(p->H[0],p->N[0],p->CA[0],p->C[0])-M_PI;
	int n=p->n_res-1;
	p->psi[n]=calc_dihedralf_angle(p->N[n],p->CA[n],p->C[n],p->O[n])-M_PI;
	for(int i=0;i<p->n_res;i++){
		for(int j=0;j<3;j++){
			b1[j]=p->N[i][j]-p->CA[i][j];
			b2[j]=p->C[i][j]-p->CA[i][j];
		}
		normalize_d(b1,3);
		normalize_d(b2,3);
		alpha[i]=acos(scal_d(b1,b2,3));
	}
	compute_dihedrals(p);
	for(int i=0;i<p->n_res;i++){
		//printf("%d alpha: %lf\n",i,alpha[i]/M_PI*180);
		error=CAT_add_peptide(p,i,  p->phi[i],M_PI-alpha[i],p->psi[i]);
	}
	compute_dihedrals(p);

	return error;
}

void print_bond_errors(FILE *stream,cat_prot *p)
{
    double
        NCa_bond[3],
        CaC_bond[3],
        CN_bond[3],
        norm_NCa,
        norm_CaC,
        norm_CN;

    fprintf(stream,"bonds test:\n");
    for(int i=0;i<p->n_res;i++)
    {
        for(int j=0;j<3;j++)
        {
            NCa_bond[j] = p->CA[i][j] -p->N [i][j];
            CaC_bond[j] = p->C [i][j] -p->CA[i][j];
            if(i<p->n_res-1)
            {
                CN_bond[j] = p->N[i+1][j] - p->C[i][j];
            }
        }
        norm_NCa = norm_d(NCa_bond, 3);
        norm_CaC = norm_d(CaC_bond, 3);
        norm_CN  = norm_d(CN_bond,  3);
        if (fabs(norm_NCa-CAT_Rbond_CaN) > BOND_MAX_DEVIATION || fabs(norm_CaC-CAT_Rbond_CCa) > BOND_MAX_DEVIATION || fabs(norm_CN-CAT_Rbond_CN) > BOND_MAX_DEVIATION)
        {
            fprintf(stream,"Nca: %g CaC: %g CN: %g\t", fabs(norm_NCa-CAT_Rbond_CaN), fabs(norm_CaC-CAT_Rbond_CCa), fabs(norm_CN-CAT_Rbond_CN));
        }
    }
    fprintf(stream,"\n");
    fflush(stream);
}

void print_joint_angles_errors(FILE *stream, cat_prot *p)
{
    double
        angle_CaCN,
        angle_CNCa;

    fprintf(stream,"joint angles test:\n");
    for(int i=1;i<p->n_res;i++)
    {
        angle_CNCa = angle_ABC( p->C [i-1], p->N[i]  , p->CA[i]);
        angle_CaCN = angle_ABC( p->CA[i-1], p->C[i-1], p->N [i]);
        if (fabs(angle_CNCa-CAT_angle_CNCa) > ANGLE_MAX_DEVIATION || fabs(angle_CaCN-CAT_angle_CaCN) > ANGLE_MAX_DEVIATION)
        {
            fprintf(stream,"%d CNCa: %g CaCN: %g \t", i,fabs(angle_CNCa-CAT_angle_CNCa), fabs(angle_CaCN-CAT_angle_CaCN));
        }
    }
    fprintf(stream,"\n");
    fflush(stream);
}
void print_omega_errors(FILE *stream, cat_prot *p)
{
    double
        omega;

    fprintf(stream,"omega test:\n");
    for(int i=1;i<p->n_res;i++)
    {
        omega = calc_dihedralf_angle(p->CA[i-1], p->C[i-1], p->N[i], p->CA[i]);
        omega = gsl_sf_angle_restrict_pos (omega);
        if (fabs(omega-M_PI) > OMEGA_MAX_DEVIATION)
        {
            fprintf(stream,"%g \n", fabs(omega-M_PI));
        }
    }
    fprintf(stream,"\n");
    fflush(stream);
}


