/**
 * @file
 * @brief Program used to evalueat efficiency of concerted rotation in constrained space of peptide backbone
 */
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sf.h>
#include <sys/types.h>
#include <unistd.h>


#include "../../lib/Caterpillar_energies.h"
#include "../../lib/my_memory.h"
#include "../../lib/my_geom.h"
#include "../../lib/geom_prop.h"
#include "../../lib/Caterpillar_IO.h"
#include "../../lib/CAT_moves.h"
#include "../../lib/quaternions.h"
#include "../../lib/histogram.h"

#define ACC 1
#define REJ 0

#define BENDING_FORCE_CONST 50000.0
#define BOND_FORCE_CONST 50000.0
#define DIHEDRAL_FORCE_CONST 50000.0

#define CATENR_CaCa_SPRING 40
// N H C O CA CB
//#define CATENR_SAW_O_O    2.7
#define CATENR_SAW_O_O    2.7
//#define CATENR_SAW_C_C    2.9
#define CATENR_SAW_C_C    0.0
//#define CATENR_SAW_CA_CA  2.8
#define CATENR_SAW_CA_CA  4.0
//#define CATENR_SAW_CB_CB  3.0
#define CATENR_SAW_CB_CB  0.0
//#define CATENR_SAW_CB_N   2.7
#define CATENR_SAW_CB_N   0.0
//#define CATENR_SAW_CB_C   2.678
#define CATENR_SAW_CB_C   0.0
//#define CATENR_SAW_CB_O   2.85
#define CATENR_SAW_CB_O   0.0
//#define CATENR_SAW_H_H    2.0
#define CATENR_SAW_H_H    2.0
//#define CATENR_SAW_O_H    2.35
#define CATENR_SAW_O_H    2.35


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

typedef struct {
    double *Hydro_T;
    double *Hydro;
    double **M;
}energy_par;


//pezzotto per rescale...
int prot_rescale ( cat_prot *p);
void rebuild ( cat_prot *p);
//


mc_traj_data * mc_traj_data_alloc(int seed,double temp);
void mc_traj_data_free( mc_traj_data * mc_trj);


double Ener_total(cat_prot *p,energy_par *ep);
void Reset_energy(mc_traj_data *mc_trj,energy_par *ep);
void Compute_delta_en(mc_move_data *mvdt,conf *C,energy_par *ep);
double Compute_energy(const cat_prot *p);
int Metropolis( mc_traj_data *mctrj );


void CAT_copy(cat_prot *dest, const cat_prot *orig);
void Init_MC(mc_move_data **mvdt, mc_traj_data **mc_traj, energy_par **ep, char file_conf[1024],char file_pot[1024]);

double protein_energy(const cat_prot *p);
double overlap_energy(const cat_prot *p);
double spring(const double force_constant, const double equil_val, const double val);
double bending_energy(const cat_prot *p);
double bond_energy(const cat_prot *p);
double dihedral_energy(const cat_prot *p);
int metropolis(const double beta, const double E, const double E_new, gsl_rng *rng_r);

#define NUM_OF_RESIDUES 20

int main(int argc, char *argv[])
{
    int
        maxIter=10e7,
        k=0,
        concerted=0,
        accepted=0,
        rejected=0;

    double
        sigma=0.05,
        beta=1.0,
        *point;
    point = (double*)calloc(2, sizeof(double));

    char
        histogramName[1024];

    FILE
        *errFile=fopen("xerr.dat","w"),
        *dihed_out=fopen("dihedrals.dat","w");

    energy_par
        *en_par;

    mc_traj_data
        *mc_trj;

//------------------------------------------
    double
        E,
        test_E,
        *origin,
        *dihed;

    origin = (double*)calloc(3, sizeof(double));
    for(int i=0;i<3;i++){origin[i]=100.0;}

    dihed = (double*)calloc(NUM_OF_RESIDUES*2, sizeof(double));
    for(int i=0;i<NUM_OF_RESIDUES;i++){dihed[2*i]=M_PI;dihed[2*i+1]=M_PI;}

    char
        *seq;

    seq = (char*)calloc(NUM_OF_RESIDUES, sizeof(char));
    for(int i=0;i<NUM_OF_RESIDUES;i++){seq[i]='A';}

    gsl_rng
        *rng_r;

    rng_r=gsl_rng_alloc(gsl_rng_taus2);
    gsl_rng_set(rng_r,1216);

    mc_move_data
        *mv;
    mv = CATMV_mc_move_data_alloc(NUM_OF_RESIDUES,NUM_OF_RESIDUES);

    printf("Doing dihedral sampling.\n");

    cat_prot 
        *protein = NULL,
        *test_protein = NULL;
    protein      = CAT_build_from_dihed( NUM_OF_RESIDUES, 5, origin, dihed, seq);
    test_protein = CAT_build_from_dihed( NUM_OF_RESIDUES, 5, origin, dihed, seq);

//------------------------------------------


    histogram
        *psi_phi         = histogram_init  (-M_PI, M_PI+0.000001, 180, 0.1);
    histogram_add_dimension(-M_PI, M_PI+0.000001, 180, psi_phi );

    mc_move_data
        *mc_mvdt;
    Init_MC(&mc_mvdt, &mc_trj, &en_par, "unused",argv[1]);

    // Initialize energy
    E=test_E=protein_energy(protein);

    //MAIN LOOP
    for(int i=0;i<maxIter;i++)
    {
        // geometry check
        if(i>0 && i%1000==0)
        {
            fprintf(errFile, "%d --\n",i);
            print_bond_errors(errFile, protein, __FILE__, __LINE__);
            print_joint_angles_errors(errFile, protein, __FILE__, __LINE__);
            print_omega_errors(errFile, protein, __FILE__, __LINE__);
        }

        // Histogram sampling
        if(i>10000 && i%50==0)
        {
            for (size_t resi=1; resi < protein->n_res-1; resi++)
            {
                point[0]=protein->psi[resi];
                point[1]=protein->phi[resi];
                histogram_add(point, psi_phi);
                fprintf(dihed_out, "%lf %lf  ", protein->phi[resi], protein->psi[resi]);
            }
            fprintf(dihed_out,"\n");
        }

        // Print out histogram and protein in pdb
        if(i>5000 && i%10000==0)
        {
            snprintf(histogramName, sizeof(histogramName), "data_%03.12lf_time_%08i.dat", sigma, i);
            histogram_print(histogramName, "matrix2d", psi_phi);

            CATIO_cat2pdb("prova.pdb", "a", "--", protein, 1);
            printf("\n\nAcceptedMoves: %d RejectedMoves: %d\n\n", accepted, rejected);
        }

        // MOVES
        if(i<5000) //first 5000 steps equil just with pivot move
        {
            CATMV_pivot(mv, test_protein, rng_r, protein->n_res/2);
        } else if (gsl_rng_uniform(mc_trj->rng_r)>0.5) {
            CATMV_pivot(mv, test_protein, rng_r, protein->n_res/5);
        } else {
            CATMV_concerted_rot(mv, test_protein, rng_r, sigma);
            concerted=1;
        }


        test_E=protein_energy(test_protein);

        if (metropolis(beta, E, test_E, rng_r))
        {
            // Move accepted
            CAT_copy(protein, test_protein);
            E=test_E;
            if (concerted){
                accepted += 1;
            }
        }
        else 
        {
            // Move rejected
            CAT_copy(test_protein, protein);
            test_E=E;
            if (concerted){
                rejected += 1;
            }
        }
        concerted=0;

    } // END MAIN LOOP

    printf("\n\nAcceptedMoves: %d RejectedMoves: %d\n\n", accepted, rejected);

    // Free used memory 
    mc_traj_data_free(mc_trj);
    CATMV_mc_move_data_free(mc_mvdt);
    free(en_par->Hydro_T);
    free(en_par->Hydro);
    free(en_par->M);
    free(en_par);
    histogram_free(psi_phi);
    free(point);
    fclose(dihed_out);
    return 0;
}

double protein_energy(const cat_prot *p)
{
    double
        energy = 0.0;
    energy += overlap_energy(p);
    energy += bending_energy(p);
    energy += bond_energy(p);
    energy += dihedral_energy(p);

    return energy;
}

double overlap_energy(const cat_prot *p)
{
    double
        e_SAW=0;

    for(int i=0;i<p->n_res;i++) {
		for(int j=0;j<p->n_res;j++) {
			int d= (i>=j)? i-j:j-i;
			if(d>=2) {
                if(1){
				    e_SAW += CATENR_Saw(p->CA[i],p->CA[j], CATENR_SAW_CA_CA);
				    e_SAW += CATENR_Saw(p->CA[j],p->CA[i], CATENR_SAW_CA_CA);
				    if(e_SAW>0) {
//					    printf("A1 CA SAW clash at i=%d j=%d\n",i,j);
					    return e_SAW; 
				    }
                }
                if(1){
                    e_SAW += CATENR_Saw(p->O[i], p->H[j], CATENR_SAW_O_H);
                    e_SAW += CATENR_Saw(p->H[i], p->O[j], CATENR_SAW_O_H);
                    if(e_SAW>0) {
//                        printf("A2 O - H SAW clash at i=%d j=%d\n",i,j);
                        return e_SAW; 
                    }
                }
			} else if (i!=j) {
                if(1){
				    e_SAW += CATENR_Saw(p->O[i], p->O[j], CATENR_SAW_O_O); // O_{i} - O_{i+1, i-1}
				    e_SAW += CATENR_Saw(p->O[j], p->O[i], CATENR_SAW_O_O); // O_{i} - O_{i+1, i-1}
				    if(e_SAW>0) {
//					    printf("A3 O - O SAW clash at i=%d j=%d\n",i,j);
					    return e_SAW; 
				    }
                }
                if(1){
				    e_SAW += CATENR_Saw(p->H[i], p->H[j], CATENR_SAW_H_H); // H_{i} - H_{i+1, i-1}
				    e_SAW += CATENR_Saw(p->H[j], p->H[i], CATENR_SAW_H_H); // H_{i} - H_{i+1, i-1}
				    if(e_SAW>0) {
//					    printf("A4 H - H SAW clash at i=%d j=%d\n",i,j);
					    return e_SAW; 
				    }
                }
                if(1){
				    e_SAW += CATENR_Saw(p->O[i], p->H[j], CATENR_SAW_O_H);
				    e_SAW += CATENR_Saw(p->H[i], p->O[j], CATENR_SAW_O_H);
				    if(e_SAW>0) {
//					    printf("A5 O - H SAW clash at i=%d j=%d\n",i,j);
					    return e_SAW; 
				    }
                }
			} else if (i==j) {
                if(1){
				    e_SAW += CATENR_Saw(p->O[i], p->H[j], CATENR_SAW_O_H);
				    e_SAW += CATENR_Saw(p->O[j], p->H[i], CATENR_SAW_O_H);
				    if(e_SAW>0) {
//					    printf("A6 O - H SAW clash at i=%d j=%d\n",i,j);
					    return e_SAW; 
				    }
                }
			}
		}
	}
    return e_SAW;
}

double spring(const double force_constant, const double equil_val, const double val)
{
    return force_constant * (equil_val-val) * (equil_val-val);
}

double bending_energy(const cat_prot *p)
{
    double
        bending_E = 0.0;

	for(int i = 0; i < p->n_res; i++)
    {
        bending_E += spring(BENDING_FORCE_CONST, CAT_angle_NCaC, angle_ABC(p->N[i] , p->CA[i], p->C[i]));
        bending_E += spring(BENDING_FORCE_CONST, CAT_angle_CaCO, angle_ABC(p->CA[i], p->C[i] , p->O[i]));
        if(i < (p->n_res - 1))
        {
            bending_E += spring(BENDING_FORCE_CONST, CAT_angle_CaCN, angle_ABC(p->CA[i], p->C[i]  , p->N[i+1]));
            bending_E += spring(BENDING_FORCE_CONST, CAT_angle_CNH , angle_ABC(p->C[i] , p->N[i+1], p->H[i+1]));
        }
        if(i > 0)
        {
            bending_E += spring(BENDING_FORCE_CONST, CAT_angle_CNCa, angle_ABC(p->C[i-1], p->N[i], p->CA[i]));
        }
    }
    return bending_E;
}

double bond_energy(const cat_prot *p)
{
    double
        bond_E = 0.0;

	for(int i = 0; i < p->n_res; i++)
    {
        bond_E +=  spring(BOND_FORCE_CONST, CAT_Rbond_NH , dist_d(p->N[i] , p->H[i] , 3));
        bond_E +=  spring(BOND_FORCE_CONST, CAT_Rbond_CaN, dist_d(p->CA[i], p->N[i] , 3));
        bond_E +=  spring(BOND_FORCE_CONST, CAT_Rbond_CCa, dist_d(p->C[i] , p->CA[i], 3));
        bond_E +=  spring(BOND_FORCE_CONST, CAT_Rbond_CO , dist_d(p->C[i] , p->O[i] , 3));
        if(i < (p->n_res - 1))
        {
            bond_E += spring(BOND_FORCE_CONST, CAT_Rbond_CN , dist_d(p->C[i] , p->N[i+1] , 3));
        }
    }
    return bond_E;
}

double dihedral_energy(const cat_prot *p)
{
    double
        dihedral_E = 0.0;

	for(int i = 0; i < p->n_res; i++)
    {
        if (i < (p->n_res - 1))
        {
            dihedral_E += spring(DIHEDRAL_FORCE_CONST, M_PI , dihedralangle_ABCD(p->O[i], p->C[i], p->N[i+1], p->H[i+1]));
        }
    }
    return dihedral_E;
}

double Ener_total(cat_prot *p,energy_par *ep)
{
	double e_SAW=0;
	double energ_CaCa=0;
	double energ_HB=0;
	double energ_Wat=0;
	double energ_Bend=0;
	for(int i=0;i<p->n_res;i++) {
		for(int j=0;j<p->n_res;j++) {
			int d= (i>=j)? i-j:j-i;
			if(d>2) {
				e_SAW += CATENR_Saw(p->CA[i],p->CA[j], CATENR_SAW_CA_CA);
				e_SAW += CATENR_Saw(p->O[i], p->H[j], CATENR_SAW_O_H);
				e_SAW += CATENR_Saw(p->H[i], p->O[j], CATENR_SAW_O_H);
				if(e_SAW>0) {
					printf("AHHH CA SAW clash at i=%d j=%d\n",i,j);
					return e_SAW; 
				}
			} else if (i!=j) {
				e_SAW += CATENR_Saw(p->O[i], p->O[j], CATENR_SAW_O_O); // O_{i} - O_{i+1, i-1}
				e_SAW += CATENR_Saw(p->H[i], p->H[j], CATENR_SAW_H_H); // H_{i} - H_{i+1, i-1}

				e_SAW += CATENR_Saw(p->O[i], p->H[j], CATENR_SAW_O_H);
				e_SAW += CATENR_Saw(p->H[i], p->O[j], CATENR_SAW_O_H);
			} else if (i==j && p->CB!= NULL) {
				e_SAW += CATENR_Saw(p->O[i], p->H[j], CATENR_SAW_O_H);
			}
			if(e_SAW>0) {
				printf("AHHH OTHER SAW clash at i=%d j=%d\n",i,j);
				return e_SAW; 
			}
		}
		energ_Wat	+=CATENR_Wat ( p,i,CAT_S-1,ep->Hydro,ep->Hydro_T,CATENR_WAT_PREF);
		energ_Bend+=CATENR_Ca_Ca_Bend(p->C[i],p->N[i],CAT_Rbond_NC,CATENR_CaCa_SPRING);
	}
	//return 0.5*energ_CaCa+0.5*energ_HB+energ_Wat+energ_Bend;
	return energ_Bend+e_SAW;
}


mc_traj_data * mc_traj_data_alloc(int seed,double temp)
{
	mc_traj_data *mctrj=(mc_traj_data*)malloc(sizeof(mc_traj_data));
	mctrj->step=0;
	mctrj->rng_r=gsl_rng_alloc(gsl_rng_taus2);
	gsl_rng_set(mctrj->rng_r,seed);
	mctrj->beta=1/temp; //Kb=1
	mctrj->temp=temp;
	mctrj->old.p=NULL;
	mctrj->old.E=-1;
	mctrj->new.p=NULL;
	mctrj->new.E=-1;
	return mctrj;
}

void mc_traj_data_free( mc_traj_data * mctrj)
{
	if(mctrj!=NULL)
	{
		gsl_rng_free(mctrj->rng_r);
        CAT_prot_free(mctrj->old.p);
        free(mctrj->old.Dcontacts);
        CAT_prot_free(mctrj->new.p);
        free(mctrj->new.Dcontacts);
		free(mctrj);
	}
}


void CAT_copy(cat_prot *dest, const cat_prot *orig)
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
				"CAT_copy. Proteins have different numbers of atoms. dest: %d. orig %d",
				dest->n_atoms,orig->n_atoms);
		failed(err);
	}
	if(dest->n_atom_per_res != orig->n_atom_per_res)
	{
		sprintf(err,
				"CAT_copy. Proteins have different numbers of atoms per residue. dest: %d. orig %d",
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


void Init_MC(mc_move_data **mvdt, mc_traj_data **mc_traj, energy_par **ep, char file_conf[1024],char file_pot[1024])
{
	//int N_res=4; //Length of the protein
	int N_res=20; //Length of the protein
	int N_atoms_per_res=5;
	int seed=getpid();
	double temp=2.0;
	double C[3]={0.,0.,0.};
	mc_traj_data *mct;
	mc_move_data *mv=CATMV_mc_move_data_alloc(N_res,N_res);
	energy_par * en = (energy_par*)malloc(sizeof(energy_par));
	mct=mc_traj_data_alloc(seed,temp);
	mct->old.p=CAT_prot_alloc(N_res,N_atoms_per_res);  //n_atom=5: no CB. put n_atom=6 to include them (untested)
	CAT_set_prot_linear(mct->old.p,C,0.0);
	CAT_set_residues_fasta(mct->old.p,20,"AAAAAAAAAAAAAAAAAAAA");
	//CAT_set_residues_fasta(mct->old.p,4,"AAAA");
//  CAT_set_residues_fasta(mct->old.p,9,"AAAAAAAAA");
	//CAT_set_residues_fasta(mct->old.p,80,"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"); //input sequence..
	//I have still to implement reading a configuration from file. It is quite 
	//easy but it is not necessary for the moment.
	//CATIO_pdb2cat (&mc_traj.old.p,&n_prot,5,ATnames_to_read,'A',argv[1]);
	//
	//alloc new protein, copy old over it.
	mct->new.p=CAT_prot_alloc(N_res,N_atoms_per_res);
	CAT_copy(mct->new.p,mct->old.p);
	//energy parameters and energy computation
	double *HT;
	double *H;
	double **M;
	CATENR_read_energy_params(CAT_S-1,&H,&HT,&M,file_pot);
	en->Hydro_T	=HT;
	en->Hydro 	=H;
	en->M     	=M;
	mct->old.E=Ener_total(mct->old.p,en);
	mct->old.Dcontacts=d1t(mct->old.p->n_res);
	mct->old.Ewat=CATENR_Wat_total ( mct->old.p,en->Hydro,en->Hydro_T,CATENR_WAT_PREF);

	mct->new.E=mct->old.E;
	mct->new.Dcontacts=d1t(mct->new.p->n_res);
	mct->new.Ewat=mct->old.Ewat;
	//
	*ep=en;
	*mvdt=mv;
	*mc_traj=mct;
}

void Reset_energy(mc_traj_data *mc_trj,energy_par *ep)
{
	double E=Ener_total(mc_trj->old.p,ep);
	fprintf(stderr,"--Energy check:st %d E= %16.10lf; En= %20.14lf; D = %20.14lf\n",mc_trj->step,E,mc_trj->new.E,E-mc_trj->new.E);
	//CATIO_cat2pdb("tua_mamma_vacca.pdb","a","--",mc_trj->old.p,1);
	//prot_rescale(mc_trj->old.p);
	//CATIO_cat2pdb("tua_mamma_vacca.pdb","a","--",mc_trj->old.p,1);
	//E=Ener_total(mc_trj->old.p,ep);
	//prot_rescale(mc_trj->new.p);
	mc_trj->old.E=E;
	mc_trj->old.Ewat=CATENR_Wat_total ( mc_trj->old.p,ep->Hydro,ep->Hydro_T,CATENR_WAT_PREF);
	//mc_trj->old.Ewat=0.0;
	mc_trj->new.E=Ener_total(mc_trj->new.p,ep);
	mc_trj->new.Ewat=mc_trj->old.Ewat;
}

void Compute_delta_en(mc_move_data *mvdt,conf *C, energy_par *ep)
{
	int n=mvdt->N_moved;
	int N=mvdt->N_pairs;
	int r1=-1;
	int r2=-1;
	cat_prot *p= C->p;
	double e_SAW=0;
	double energ_CaCa=0;
	double energ_HB=0;
	double energ_Wat=0;
	double energ_Bend=0;
	double r,bond;
	//bending energy only for termini (I am computing it also for pivots..
	//this needs to be corrected in the future
	if(0==n)
	{
		C->DE=0;
		return;
	}
	r1=mvdt->moved_res[0];
	r2=mvdt->moved_res[n-1];
	energ_Bend+=CATENR_Ca_Ca_Bend(p->C[r1],p->N[r1],CAT_Rbond_NC,CATENR_CaCa_SPRING);
	energ_Bend+=CATENR_Ca_Ca_Bend(p->C[r2],p->N[r2],CAT_Rbond_NC,CATENR_CaCa_SPRING);
	//loop on the interacting atom pairs (needs to be optimized!!)
	int step=1;
	for(int i=0;i<p->n_res;i++) {
		C->Dcontacts[i]=0;
	}
	for(int i=0;i<N;i++) {
		r1=mvdt->mod_pairs[i].i_1;
		r2=mvdt->mod_pairs[i].i_2;
		int d=r1>=r2? r1-r2:r2-r1;
		if(d>2) {
			r =dist_d(p->CA[r1],p->CA[r2],3);
			e_SAW  = CATENR_Saw_eq (r, CATENR_SAW_CA_CA);
			e_SAW += CATENR_Saw		 (p->O[r1], p->H[r2], CATENR_SAW_O_H);
			e_SAW += CATENR_Saw		 (p->H[r1], p->O[r2], CATENR_SAW_O_H);

			if(e_SAW>0) {
			 	C->DE=energ_CaCa+energ_HB+energ_Bend+e_SAW;
				return;
			}
			else if(r<CATENR_CaCa_CUTOFF) {//This is the largest cut-off, so we keep it as a reference
				bond=CATENR_bond(r);
				energ_CaCa+=CATENR_Ca_Ca_eq(p->residues[r1],p->residues[r2],bond,ep->M,CATENR_CaCa_PREF);
				//necessary to compute the water exposure..
				C->Dcontacts[r1]+=bond;
				C->Dcontacts[r2]+=bond;
				//this is possibly suboptimal since the cutoff is much lower. The same 
				//goes for the SAW, but for the moment who cares.
				energ_HB	+=CATENR_HB 	(p->H[r1],p->N[r1],p->O[r2],p->C[r2],CATENR_HB_PREF);
				energ_HB	+=CATENR_HB 	(p->H[r2],p->N[r2],p->O[r1],p->C[r1],CATENR_HB_PREF);
			}
		} else if (r1 != r2) {
			e_SAW += CATENR_Saw(p->O[r1], p->O[r2], CATENR_SAW_O_O); // O_{i} - O_{i+1, i-1}
			e_SAW += CATENR_Saw(p->H[r1], p->H[r2], CATENR_SAW_H_H); // H_{i} - H_{i+1, i-1}

			e_SAW += CATENR_Saw(p->O[r1], p->H[r2], CATENR_SAW_O_H);
			e_SAW += CATENR_Saw(p->H[r1], p->O[r2], CATENR_SAW_O_H);
		} else if(r1==r2) {
			e_SAW += CATENR_Saw(p->O[r1], p->H[r2], CATENR_SAW_O_H);
		}
		if(e_SAW>0) {
			break;
		}
	}
	//C->DE=energ_CaCa+energ_HB+energ_Bend+e_SAW;
 	C->DE=energ_Bend+e_SAW;
}


void Compute_energy_new(mc_move_data *mvdt,mc_traj_data *mctrj,energy_par *ep)
{
	double ener_w;
	double Ewat=0;
	int type;
	Compute_delta_en(mvdt,&mctrj->old,ep);
	Compute_delta_en(mvdt,&mctrj->new,ep);
	//water exposure...
	for(int i=0;i<mctrj->old.p->n_res;i++)
	{
		mctrj->new.p->contacts[i]=mctrj->old.p->contacts[i]-mctrj->old.Dcontacts[i]+mctrj->new.Dcontacts[i];
		type=mctrj->old.p->residues[i];
		ener_w=ep->Hydro[type]*(ep->Hydro_T[type]-mctrj->new.p->contacts[i]);
		ener_w=ener_w>0?ener_w:0;
		Ewat+=ener_w;
	}
	//Ewat*=CATENR_WAT_PREF;
	mctrj->new.Ewat=CATENR_WAT_PREF*Ewat; //dannato prefisso!
	mctrj->new.E=mctrj->old.E-mctrj->old.DE + mctrj->new.DE - mctrj->old.Ewat + mctrj->new.Ewat;
	//mctrj->new.E=mctrj->old.E-mctrj->old.DE + mctrj->new.DE;
}

int metropolis(const double beta, const double E, const double E_new, gsl_rng *rng_r)
{
    double randomVal = gsl_rng_uniform(rng_r);
//    printf("E: %g | test_E: %g | beta: %g | log(rand): %g | beta*(E-E_new)=%g\n", E, E_new, beta, log(randomVal), (beta * (E-E_new) ));
    if ((beta * (E-E_new)) > log(randomVal)){
//        printf("ACC\n");
        return ACC;
    }else{
//        printf("REJ\n");
        return REJ;
    }
//    return ((beta * (E-E_new) ) > log(gsl_rng_uniform(rng_r))) ? ACC : REJ;
}

int Metropolis( mc_traj_data *mctrj )
{
	mctrj->step+=1;
	double DE=mctrj->old.E-mctrj->new.E;
//	double DE=Compute_energy(mctrj->old.p)-Compute_energy(mctrj->new.p);
	double LR=log(gsl_rng_uniform(mctrj->rng_r));
	//fprintf(stderr,"--Metro: DE = %lf LR = %lf\n",DE,LR);
	if(mctrj->beta*DE>LR)
	{
		CAT_copy(mctrj->old.p,mctrj->new.p);
		mctrj->old.E=mctrj->new.E;
		mctrj->old.Ewat=mctrj->new.Ewat;
		return ACC;
	}
	else
	{
		CAT_copy(mctrj->new.p,mctrj->old.p);
		mctrj->new.E=mctrj->old.E;
		mctrj->new.Ewat=mctrj->old.Ewat;
		return REJ;
	}
}

int prot_rescale ( cat_prot *p)
{
	int error;
	double r;
	double alpha[p->n_res];
	double b1[3],b2[3];
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
	p->phi[0]=dihedralangle_ABCD(p->H[0],p->N[0],p->CA[0],p->C[0])-M_PI;
	int n=p->n_res-1;
	p->psi[n]=dihedralangle_ABCD(p->N[n],p->CA[n],p->C[n],p->O[n])-M_PI;
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

void rebuild(cat_prot *p)
{
	for(int i = 0; i < p->n_res; i++)
	{
		CAT_add_peptide(p, i, p->phi[i],M_PI-CAT_angle_NCaC,p->psi[i]);
	}
	return;
}


