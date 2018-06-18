/**
 * @file
 * @brief Minimal implementation of concerted rotation move 
 * Main purpose is to enable easier implementation in any MC code.
 */

#include <stdio.h>
#include <math.h>
#include <unistd.h>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>

#include <sys/types.h>
#include "generic_move/CR_precomp.h"

#define MOVE_LENGHTH 7



/*******PARAMETERS OF THE MODEL*******/
//Backbone bond lengths
#define CAT_Rbond_CaN 1.4500000 //sist->Rbond[0]= 1.45000;
#define CAT_Rbond_CCa 1.5200000	//sist->Rbond[1]=1.52000;
#define CAT_Rbond_CN  1.3300000	//sist->Rbond[3]= 1.33000;
#define CAT_Rbond_NC  2.4479800	//sist->Rbond[5], used for Bending instead of the angle..

//Backbone angles
#define CAT_angle_CaCN 2.017600615305445
#define CAT_angle_CNCa 2.1275563581810877
#define CAT_angle_NCaC 1.9373154697137058
/*^^^^^^^PARAMETERS OF THE MODEL^^^^^^^*/

// structure used in solver and T7_solver functions
struct rparams {
	double s;
	gsl_matrix *B;
	gsl_vector *wspace;
	cr_input_data bb_in;
	cr_input_data bb_out;
};

// functions needed for concerted rotation move (have to be included in every implementation!)
int random_rot(cr_input_data bb_out, cr_input_data bb_in, gsl_rng *rng_r, double sigma);
int solver(cr_input_data bb_out, cr_input_data bb_in, int angle,double delta_mov);
double jac_Det(int a, cr_input_data bb_in);
int T7_solver (const gsl_vector *x, void *p, gsl_vector *f);
gsl_matrix * gram_schmidt ( gsl_vector *v);


// DO NOT FORGET TO INCLUDE ALSO message.[c,h] due to usage in CR_precomp.[c,h] !!!

gsl_rng *rng_r;
int main(int argc, char *argv[])
{
	int
        N_res = 20,             // number of "residues"
        c,                      // number that define the first "residue" from which concerted rotation is done
        maxIteration = 1e8,     // how many concerted rotation are done
        error;                  // value to catch error code from random_rot

    double
        dihedrals[N_res][2],    // [x][0] is 
        sigma = 0.10;           // fixed value of sigma used in concerted rotation move (see Documentation)


    FILE
        *output;
        output = fopen ("dihedrals.dat", "w+");

    // alloc interface for random rotation
    cr_input_data // see CR_precomp.h
        bb_in,
        bb_out;

    alloc_cr_input_data(&bb_in);
    alloc_cr_input_data(&bb_out);

    // Backbone joint angles and bond lengths are fixed and defined here
    // Check the Mathematica notebook for details
    //
    // If more then 3 residues are supposed to be moved then more similar
    // lines have to be included.
    gsl_vector_set(bb_in.r,0,CAT_Rbond_CaN); 
    gsl_vector_set(bb_in.r,1,CAT_Rbond_CCa-CAT_Rbond_NC*cos(CAT_angle_CaCN)); 
    gsl_vector_set(bb_in.r,2,CAT_Rbond_CaN);
    gsl_vector_set(bb_in.r,3,CAT_Rbond_CCa-CAT_Rbond_NC*cos(CAT_angle_CaCN)); 
    gsl_vector_set(bb_in.r,4,CAT_Rbond_CaN);
    gsl_vector_set(bb_in.r,5,CAT_Rbond_CCa-CAT_Rbond_NC*cos(CAT_angle_CaCN)); 
    gsl_vector_set(bb_in.r,6,CAT_Rbond_CaN);

    gsl_vector_set(bb_in.d,0,0);
    gsl_vector_set(bb_in.d,1,CAT_Rbond_CN*sin(CAT_angle_CaCN));
    gsl_vector_set(bb_in.d,2,0);
    gsl_vector_set(bb_in.d,3,CAT_Rbond_CN*sin(CAT_angle_CaCN));
    gsl_vector_set(bb_in.d,4,0);
    gsl_vector_set(bb_in.d,5,CAT_Rbond_CN*sin(CAT_angle_CaCN));
    gsl_vector_set(bb_in.d,6,0);

    gsl_vector_set(bb_in.bend_angles,0,CAT_angle_NCaC);
    gsl_vector_set(bb_in.bend_angles,1,CAT_angle_CNCa-CAT_angle_CaCN);
    gsl_vector_set(bb_in.bend_angles,2,CAT_angle_NCaC);
    gsl_vector_set(bb_in.bend_angles,3,CAT_angle_CNCa-CAT_angle_CaCN);
    gsl_vector_set(bb_in.bend_angles,4,CAT_angle_NCaC);
    gsl_vector_set(bb_in.bend_angles,5,CAT_angle_CNCa-CAT_angle_CaCN);
    gsl_vector_set(bb_in.bend_angles,6,CAT_angle_NCaC);

    // copy bb_in to bb_out (see CR_precomp.[h,c])
    memcpy_cr_input_data(&bb_out, &bb_in);

    //random number generator initialization from first argument of the program
    rng_r = gsl_rng_alloc(gsl_rng_taus2);
    gsl_rng_set(rng_r, atoi(argv[1]));


    //create random dihedral configuration
    for(int i = 0; i < N_res; i++)
    {
        dihedrals[i][0]=2*M_PI*gsl_rng_uniform(rng_r)-M_PI;
        dihedrals[i][1]=2*M_PI*gsl_rng_uniform(rng_r)-M_PI;
    }

    for(int k = 0; k < maxIteration; k++) // simulation loop at given sigma
    {
        // select random part of chain
        c = 1 + gsl_rng_uniform_int(rng_r, N_res - 2 );

        // Print out statistics
        if(k%1000==0)
        {
            printf("Iteration: %i\n", k);fflush(stdout);
            for(int i = 0; i < N_res; i++)
            {
                fprintf(output, "%f\t", dihedrals[i][0]);
                fprintf(output, "%f\t", dihedrals[i][1]);
            }
            fprintf(output, "\n");
        }

        //params - PHI matrix
        for(int i=0;i<4;i++)
        {
            int j=i*2;
            //torsion angles (dihedrals)
            gsl_vector_set(bb_in.dihed_angles, j, dihedrals[c+i][0]);
        }

        //params - PSI matrix
        for(int i=0;i<3;i++)
        {
            int j=2*i+1;
            //torsion angles (dihedrals)
            gsl_vector_set(bb_in.dihed_angles,j,dihedrals[c+i][1]);
        }

        gsl_vector_memcpy ( bb_out.dihed_angles, bb_in.dihed_angles);

        // Actual move concerted rotation move
        error = random_rot(bb_out, bb_in,rng_r, sigma);
        if (error != 0){continue;}

        // Update the dihedrals after move
        dihedrals[c][0]  =gsl_sf_angle_restrict_symm(gsl_vector_get(bb_out.dihed_angles,0));
        dihedrals[c][1]  =gsl_sf_angle_restrict_symm(gsl_vector_get(bb_out.dihed_angles,1));
        dihedrals[c+1][0]=gsl_sf_angle_restrict_symm(gsl_vector_get(bb_out.dihed_angles,2));
        dihedrals[c+1][1]=gsl_sf_angle_restrict_symm(gsl_vector_get(bb_out.dihed_angles,3));
        dihedrals[c+2][0]=gsl_sf_angle_restrict_symm(gsl_vector_get(bb_out.dihed_angles,4));
        dihedrals[c+2][1]=gsl_sf_angle_restrict_symm(gsl_vector_get(bb_out.dihed_angles,5));
        dihedrals[c+3][0]=gsl_sf_angle_restrict_symm(gsl_vector_get(bb_out.dihed_angles,6));
    }

    fclose(output);
    return error;
}

int random_rot(cr_input_data bb_out, cr_input_data bb_in, gsl_rng *rng_r, double sigma)
{
	int
        error = 1,
        a,
        i,
        m=MOVE_LENGHTH,
        angles_idx[MOVE_LENGHTH];

    for (int ii = 0; ii < MOVE_LENGHTH; ii++){
        angles_idx[ii]=ii;
    }

	double
        ds=gsl_ran_gaussian(rng_r,sigma);
	//Try all 7 rotation angles.
	do{
		i=gsl_rng_uniform_int (rng_r, m);
		a=angles_idx[i];
		error=solver(bb_out,bb_in,a,ds);
		if(error)
		{
			m--;
			angles_idx[i]=angles_idx[m];
		}
	}while(m && error!=GSL_SUCCESS);
	//Work out the correct acceptance
	double ds1;
	double xi_new_ar[MOVE_LENGHTH];
	double delta_xi_ar[MOVE_LENGHTH];
	gsl_vector_view T_xi_new_v =gsl_vector_view_array(xi_new_ar,MOVE_LENGHTH);
	gsl_vector *T_xi_new=&T_xi_new_v.vector;
	gsl_vector_view delta_xi_v =gsl_vector_view_array(delta_xi_ar,MOVE_LENGHTH);
	gsl_vector *delta_xi=&delta_xi_v.vector;
	//To be improved..
	if(error==GSL_SUCCESS)
	{
		double Jo,Jn;
		Jo=fabs(jac_Det(a,bb_in));
		m=MOVE_LENGHTH;
		int status=1;
		int angles_idx_B[MOVE_LENGHTH];
        for (int ii=0;ii<MOVE_LENGHTH;ii++){angles_idx_B[ii]=i;}
		do{
			i=gsl_rng_uniform_int (rng_r, m);
			a=angles_idx_B[i];
            status=(*TmT[a])(T_xi_new,bb_out);
			if(status)
			{
				m--;
				angles_idx_B[i]=angles_idx_B[m];
			}
		}while(m && status);
		Jn=fabs(jac_Det(a,bb_out));
		if(status!=0) {
			printf("phi db -Jacobian not invertible. dice: %d Error code %d\n",a,status);
			//return status;
			error = GSL_FAILURE;
		}
		//step 2. project the difference xi_old-xi_new
		if(error==GSL_SUCCESS)
		{
			double w1,w2;
			for(i=0;i<MOVE_LENGHTH;i++)
			{
				w1=gsl_vector_get(bb_out.dihed_angles,i);
				w2=gsl_vector_get(bb_in.dihed_angles,i);
				w1=w2-w1;
				w1=GSL_MIN_DBL(w1,2*M_PI-w1);
				gsl_vector_set(delta_xi,i,w1);
			}
			gsl_blas_ddot(T_xi_new,delta_xi,&ds1);
			double norm=gsl_blas_dnrm2 (T_xi_new);
			ds1 =ds1/norm;
			//step 3. compute acceptance
			double acc= exp((ds*ds-ds1*ds1)/(2*sigma*sigma))*Jo/Jn;
			if(acc <1 && gsl_rng_uniform(rng_r)>=acc) //reject the move.
			{
				//printf("REJphi ds1 %lf ds %lf diff %lf\n",ds1*ds1,ds*ds,ds*ds-ds1*ds1);
				error = -100; 
			}
		}
	}

	if(error!=GSL_SUCCESS)
	{
		gsl_vector_memcpy(bb_out.dihed_angles,bb_in.dihed_angles);
	}
	return error;
}


int solver(cr_input_data bb_out, cr_input_data bb_in, int angle, double delta_mov)
{
    double
        T_ar[MOVE_LENGHTH];

    gsl_vector_view
        Tv = gsl_vector_view_array(T_ar, MOVE_LENGHTH);

    gsl_vector
        *T = &Tv.vector;

    gsl_matrix
        *B;

    //solver
    int
        status;

    size_t
        iter = 0;

    const size_t
        n = MOVE_LENGHTH-1;

    const gsl_multiroot_fsolver_type
        *Ts;

    gsl_multiroot_fsolver
        *s;

    double
        x_ar[MOVE_LENGHTH-1];

    gsl_vector_view
        x_v = gsl_vector_view_array(x_ar,MOVE_LENGHTH-1);

    gsl_vector
        *x = &x_v.vector;

    //function parameters
    struct rparams
        p;

    //choose an angle, get tangential vector to manifold
    status = (*TmT[angle])(T, bb_in); // for TmT see CR_precomp.[h,c]

    if(status!=0) {
        //printf("Jacobian not invertible. dice: %d Error code %d\n",angle,status);
        return GSL_FAILURE; //again, I know it is not gsl's fault..
    }
    // get orthonormalbasis to T manifold
    B = gram_schmidt (T);

    //parameters
    p.s = delta_mov;
    p.B = B;

    double
        ws_ar[MOVE_LENGHTH];

    gsl_vector_view
        ws_v = gsl_vector_view_array(ws_ar,MOVE_LENGHTH);

    p.wspace = &ws_v.vector;
    //
    p.bb_out=bb_out;
    p.bb_in	=bb_in;
    // !!!CHECK MOVE_LENGHTH-1!!!
    gsl_multiroot_function
        f = {&T7_solver, n, &p};

    gsl_vector_set_zero(x);

    //solver inizialization
    Ts = gsl_multiroot_fsolver_hybrid;
    s = gsl_multiroot_fsolver_alloc (Ts, MOVE_LENGHTH-1);
    gsl_multiroot_fsolver_set (s, &f, x);

    do
    {
        iter++;
        status = gsl_multiroot_fsolver_iterate (s);
        if (status)   // check if solver is stuck
        {
            break;
        }
        status = gsl_multiroot_test_residual (s->f, 1e-12);
    } while (status == GSL_CONTINUE && iter < 100);

    gsl_multiroot_fsolver_free (s);
    gsl_matrix_free (B);
    return status;
}

double jac_Det(int a, cr_input_data bb_in)
{
    int
        sign;

    double
        J;

    gsl_matrix
        *M;
    gsl_permutation
        *p = gsl_permutation_alloc (MOVE_LENGHTH-1);

    M = (*jac[a])(bb_in);
	gsl_linalg_LU_decomp(M,p,&sign);
	J = gsl_linalg_LU_det(M,sign);
	gsl_matrix_free(M);
	gsl_permutation_free(p);
	return J;
}

int T7_solver (const gsl_vector *x, void *p, gsl_vector *f)
{
	//get parameters
	gsl_matrix *B		=((struct rparams*)	p)->B;
	double s 				=((struct rparams*)	p)->s;
	gsl_vector *w		=((struct rparams*)	p)->wspace;
	cr_input_data bb_in =((struct rparams*) p)->bb_in;
	cr_input_data bb_out =((struct rparams*) p)->bb_out;
	gsl_vector *xi0	=bb_in.dihed_angles;
	gsl_vector *an  =bb_out.dihed_angles;

	//other vars
	int i, error;
	double w_ip1;
	//trans
	//gsl_vector_view an_=gsl_vector_subvector(an,0,7);
	gsl_vector_set(w,0,s);
	for(i=0;i<MOVE_LENGHTH-1;i++)
	{
		w_ip1=gsl_vector_get( x, i);
		gsl_vector_set(w,i+1,w_ip1);
	}
	//printf("AHHHHH scaling s =%lf\n",s);
	//multiply the coefficients s and x to obtain the angles.
	//we must set an to xi0 and sum the proposed motion to it.
	gsl_vector_memcpy(an,xi0);
	error=gsl_blas_dgemv (CblasTrans, 1., B, w, 1., an);
	//gsl_vector_set(an,6,gsl_vector_get(xi0,6));
	gsl_vector_view  w_ = gsl_vector_subvector (w, 1, MOVE_LENGHTH-1);
	T7( f, bb_out );
	T7( &w_.vector, bb_in);
	gsl_vector_sub(f,&w_.vector);
	return error;
}


gsl_matrix * gram_schmidt ( gsl_vector *v)
{
	int n=v->size;
	int i,j,k;
	double norm_v2, norm_v;
	double a;
	double vw[n][n];
	gsl_matrix *B=gsl_matrix_alloc(n,n);
	k=gsl_vector_max_index (v);
	norm_v2=0;
	//get vector
	for(i=0;i<n;i++)
	{
		vw[0][i]=gsl_vector_get(v,i);
		norm_v2+=vw[0][i]*vw[0][i];
	}
	norm_v=sqrt(norm_v2);
	//construct basis
	for(i=1;i<n;i++)
	{
		if(i==k)
		{
			vw[k][0]=1.;
			for(j=1;j<n;j++) { vw[k][j]=0; }
		}
		else
		{
			for(j=0;j<n;j++) { vw[i][j]=0.;}
			vw[i][i]=1.;
		}
	}
	//orthonormalize and save basis
	for(k=0;k<n;k++)
	{
		vw[0][k]/=norm_v;
		gsl_matrix_set(B,0,k,vw[0][k]);
	}
	for(i=0;i<n;i++)
	{
		if(i>0){
			norm_v2=0;
			for(k=0;k<n;k++) { norm_v2+=vw[i][k]*vw[i][k]; }
			norm_v=sqrt(norm_v2);
			for(k=0;k<n;k++) {
				vw[i][k]/=norm_v;
				gsl_matrix_set(B,i,k,vw[i][k]);
			}
		}
		for(j=i+1;j<n;j++)
		{
			a=0;
			for(k=0;k<n;k++) { a+=vw[j][k]*vw[i][k]; }
			for(k=0;k<n;k++) { vw[j][k]-=a*vw[i][k]; }
		}
	}
	return B;
}
