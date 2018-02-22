#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <gsl/gsl_sf.h>
#include <sys/types.h>

#include "../../lib/Caterpillar_energies.h"
#include "../../lib/my_memory.h"
#include "../../lib/my_geom.h"
#include "../../lib/geom_prop.h"
#include "../../lib/Caterpillar_IO.h"
#include "../../lib/CAT_moves.h"
#include "../../lib/quaternions.h"
#include "../../lib/histogram.h"

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


mc_traj_data * mc_traj_data_alloc(int seed,double temp);
void mc_traj_data_free( mc_traj_data * mc_trj);

gsl_rng *rng_r;
int main(int argc, char *argv[])
{
	int
        N_res=20,               // number of "residues"
        c,                      // number that define the first "residue" from which concerted rotation is done
        maxIteration=1e4,       // how many concerted rotation are done for each value of sigma
        numberOfRejects=0,      // counter for number of rejections due to move acceptance defined in article
        numberOfAllRejects=0;   // counter including also rejection due to fail in solving equtions
   
    clock_t
        begin,                  // Used to asses lenght of one concerted rotation move
        end;                    // -||-

	double
        averageTime,            // value used to store average time spend in concerted rotation move
        sigma = 0.05,           // sigma used in concerted rotation move defining sigma of normal distribution of random dispacement in Tangent space
        dihedrals[N_res][2];    // [x][0] is 

    if(argc < 2){
        printf("Have to specifie seed for random generator as argument!\n");
        return 1;
    }

    for(int i=0; i <= 50; i++)   // Loop over different sigma parameters
    {
        averageTime = 0.0;   // reset counters for each sigma
        numberOfRejects = 0; // reset
        numberOfAllRejects = 0; // reset
        sigma += 0.02 * i;   // modifie sigma

	    //alloc interface for random rotation
	    struct cr_input_data bb_in;
	    struct cr_input_data bb_out;
	    bb_in.dihed_angles=gsl_vector_alloc(7);
	    bb_in.bend_angles=gsl_vector_alloc(7);
	    bb_in.r=gsl_vector_alloc(7);
	    bb_in.d=gsl_vector_alloc(7);
	    bb_out.dihed_angles=gsl_vector_alloc(7);
	    bb_out.bend_angles=gsl_vector_alloc(7);
	    bb_out.r=gsl_vector_alloc(7);
	    bb_out.d=gsl_vector_alloc(7);

	    //backbone joint angles and bond lengths are fixed. 
	    //Check the Mathematica notebook for details.
	    for(int i=0;i<4;i++)
	    {
	    	int j=i*2;
	    	//link length
	    	gsl_vector_set(bb_in.r,j,CAT_Rbond_CaN); 
	    	//link offset
	    	gsl_vector_set(bb_in.d,i,0);
	    	//joint angles
	    	gsl_vector_set(bb_in.bend_angles,j,CAT_angle_NCaC);
	    }
	    //params - PSI matrix
	    for(int i=0;i<3;i++)
	    {
	    	int j=2*i+1;
	    	//link length
	    	gsl_vector_set(bb_in.r,j,CAT_Rbond_CCa-CAT_Rbond_NC*cos(CAT_angle_CaCN)); 
	    	//link offset
	    	gsl_vector_set(bb_in.d,i,CAT_Rbond_CN*sin(CAT_angle_CaCN));
	    	//joint angles
	    	//gsl_vector_set(cr_in.bend_angles,j,0.10995574287);
	    	gsl_vector_set(bb_in.bend_angles,j,CAT_angle_CNCa-CAT_angle_CaCN);
	    }
	    gsl_vector_memcpy ( bb_out.bend_angles, bb_in.bend_angles);
	    gsl_vector_memcpy ( bb_out.r, bb_in.r);
	    gsl_vector_memcpy ( bb_out.d, bb_in.d);


	    //random number generator
	    rng_r=gsl_rng_alloc(gsl_rng_taus2);
	    gsl_rng_set(rng_r, atoi(argv[1])); // random generator is initialized from frist argument


	    //create random dihedral configuration
	    for(int i=0;i<N_res;i++) {
	    	dihedrals[i][0]=2*M_PI*gsl_rng_uniform(rng_r)-M_PI;
	    	dihedrals[i][1]=2*M_PI*gsl_rng_uniform(rng_r)-M_PI;
	    }
	    for(int k=0;k<maxIteration;k++){ // simulation loop at given sigma
	    	c=1+gsl_rng_uniform_int (rng_r, N_res - 6 );
	    	//set input parameters for the random rotation, i.e. the torsion angles
	    	for(int i=0;i<4;i++) {
	    		int j=i*2;
	    		//torsion angles (dihedrals)
	    		gsl_vector_set(bb_in.dihed_angles,j,dihedrals[c+i][0]);
	    	}
	    	//params - PSI matrix
	    	for(int i=0;i<3;i++) {
	    		int j=2*i+1;
	    		//torsion angles (dihedrals)
	    		gsl_vector_set(bb_in.dihed_angles,j,dihedrals[c+i][1]);
	    	}
	    	gsl_vector_memcpy ( bb_out.dihed_angles, bb_in.dihed_angles);

		    //perform concerted rotation
            begin=clock();
    		int error=random_rot(bb_out, bb_in,rng_r, sigma);
            end=clock();
            averageTime+=(end-begin);
            if (error == -100){numberOfRejects++;}
            if (error != GSL_SUCCESS){numberOfAllRejects++;}

    		//update the dihedrals.
    		gsl_vector_memcpy ( bb_in.dihed_angles, bb_out.dihed_angles);
    		dihedrals[c][0]	 =gsl_sf_angle_restrict_symm(gsl_vector_get(bb_out.dihed_angles,0));
    		dihedrals[c][1]  =gsl_sf_angle_restrict_symm(gsl_vector_get(bb_out.dihed_angles,1));
    		dihedrals[c+1][0]=gsl_sf_angle_restrict_symm(gsl_vector_get(bb_out.dihed_angles,2));
    		dihedrals[c+1][1]=gsl_sf_angle_restrict_symm(gsl_vector_get(bb_out.dihed_angles,3));
    		dihedrals[c+2][0]=gsl_sf_angle_restrict_symm(gsl_vector_get(bb_out.dihed_angles,4));
    		dihedrals[c+2][1]=gsl_sf_angle_restrict_symm(gsl_vector_get(bb_out.dihed_angles,5));
    		dihedrals[c+3][0]=gsl_sf_angle_restrict_symm(gsl_vector_get(bb_out.dihed_angles,6));
            if(k%10000==0){
                printf("Iteration: %i\n", k);fflush(stdout);
            }
    	}
        printf("sigma randon_rot_time rejection_probability allRejection_probability: %lf %lf %lf %lf\n", sigma, averageTime/((double)CLOCKS_PER_SEC * maxIteration), (double)numberOfRejects/(double)maxIteration, (double)numberOfAllRejects/(double)maxIteration);
    }// end of sigma loop
	return 0;
}


