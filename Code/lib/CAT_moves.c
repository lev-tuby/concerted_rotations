/**
 * @file
 * @brief Function for manipulation with #CAT_prot strcuture
 */

#include "CAT_moves.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>
#include <stdio.h>
#include <math.h>
#include "my_memory.h"
#include "my_geom.h"
#include "../generic_move/CR_precomp.h"
//#include "Caterpillar_IO.h"

/**
 * @brief Generates D-H backbone data structure format starting from phi angle
 *
 * Parameters set starts with phi angle of first residue defined by #c value
 */
cr_input_data set_cr_params_phi (cat_prot *p, int c);

/**
 * @brief Generates D-H backbone data structure format starting from phi angle
 *
 * Parameters set starts with psi angle of first residue defined by #c value
 */
cr_input_data set_cr_params_psi(cat_prot *p, int c);

/**
 * @brief Determinant of Dini's jacobian
 * @param[in] a  			index of free variable
 * @param[in] bb_in   current point on the c.r. manifold
 */
double jac_Det (int a, cr_input_data bb_in);

/**
 * @brief Finds a point on the c.r. manifold starting from a point on its tangent space, using a root-finding algorithm.
 *
 */
int solver (cr_input_data bb_out, cr_input_data bb_in, int angle,double delta_mov);

/**
 * @brief Function passed to the root-finding algorithm.
 */
int T7_solver (const gsl_vector *x, void *p, gsl_vector *f);

/**
 * @brief Prints iteration number and state of liner solver
 */
void print_state (int iter, gsl_multiroot_fsolver * s);

/**
 * @brief Structure storing the  parameters for solving projection from tangent space to mainfold used in #T7_solver
 */
struct rparams {
	double          s;          /**< Shift along tangent space. */
	gsl_matrix      *B;         /**< Orthogonal basis set to tanget space. */
	gsl_vector      *wspace;    /**< working memory for computations */
	cr_input_data   bb_in;      /**< starting configuration in D-H convention. (4 residues) */
	cr_input_data   bb_out;     /**< new configuration in D-H convention (4 residues) */
};

/**
 *
 * #resid_pair that are moved are initialized to {-1, -1, -1.0}.
 *
 * @param[in]        max_move_size      Define how many residues could be altered in one move.
 * @param[in]        len                Length of entire protein.
 *
 * @return initialized *#mc_move_data
 */
mc_move_data * CATMV_mc_move_data_alloc(int max_move_size, int len)
{
	mc_move_data *mvdt=(mc_move_data*)malloc(sizeof(mc_move_data));
	mvdt->_l1=max_move_size;
	mvdt->_l2=len;
	int l1=mvdt->_l1;
	int l2=mvdt->_l2;
	mvdt->moved_res		=i1t(l1); //initializes to zero
	mvdt->mod_pairs=(resid_pair*)malloc(l1*l2*sizeof(resid_pair));
	resid_pair ap_default={-1,-1,-1.0};
	for(int i=0;i<l1*l2;i++) {
		mvdt->mod_pairs[i]=ap_default;
	}
	return mvdt;
}

/**
 * @param[in,out]    mvdt     #mc_move_data structure to be deallocated
 *
 * @return \c void
 */
void CATMV_mc_move_data_free( mc_move_data *mvdt)
{
	free(mvdt->moved_res);
	free(mvdt->mod_pairs);
	free(mvdt);
}

/**
 * @note The maximum cranckshaft size is 10 residues (hardcoded).
 * @note The maximum angle is 0.1 radians (hardcoded)
 *
 * @param[in,out]   *cranck_data      Move data
 * @param[in,out]   *p                protein
 * @param[in]       *rng_r            GSL random number generator state
 *
 * @return \c void
 */
void CATMV_cranck (mc_move_data *cranck_data, cat_prot *p, gsl_rng *rng_r)
{
	int k;
	int end;
	double v_1[3];

	int length=p->n_res>=10?8:p->n_res-2;
	int c=gsl_rng_uniform_int (rng_r, p->n_res-length-1);
	length=2+gsl_rng_uniform_int(rng_r,length);
	double angle= 0.2*M_PI*(-0.5+gsl_rng_uniform(rng_r));
	//check boundaries
	if(c+length >= p->n_res)
	{
		end=p->n_res-1;
		length=end-c;
	}
	else
	{
		end=c+length;
	}
	int start=c*p->n_atom_per_res+ATOM_CA+1;
	int e=end*p->n_atom_per_res+ATOM_CA; //the last CA is untouched
	int len=e-start;

	//set rotation axis
	for(k=0;k<3;k++) { v_1[k] = p->CA[end][k] -	p->CA [c][k]; }
	//rotate
	rotate(&p->coord[start],len,p->CA[c],v_1,angle);
	if(p->CB!=NULL) {
		CAT_insert_cbeta(p,end,CAT_Cb_AZIMUTH,CAT_Rbond_CCb);
		CAT_insert_cbeta(p,c,CAT_Cb_AZIMUTH,CAT_Rbond_CCb);
	}
	//compute the initial and final dihedrals too
	p->psi[c]=dihedralangle_ABCD(p->N[c],p->CA[c],p->C[c],p->N[c+1]);
	if(c>0) {
		p->phi[c]=dihedralangle_ABCD(p->C[c-1],p->N[c],p->CA[c],p->C[c]);
	}
	p->phi[end]=dihedralangle_ABCD(p->C[end-1],p->N[end],p->CA[end],p->C[end]);
	if(end<p->n_res-1) {
		p->psi[end]=dihedralangle_ABCD(p->N[end],p->CA[end],p->C[end],p->N[end+1]);
	}
	cranck_data->N_moved=length+1;
	//Save data on moved particles
	k=0;
	for(int i=0;i<cranck_data->N_moved;i++) {
		cranck_data->moved_res[i]=c+i;
		for(int j=0;j<=c;j++) {
			cranck_data->mod_pairs[k].i_1=c+i;
			cranck_data->mod_pairs[k].i_2=j;
			cranck_data->mod_pairs[k].dist=-1.0;
			k++;
		}
		for(int j=end;j<p->n_res;j++) {
			cranck_data->mod_pairs[k].i_1=c+i;
			cranck_data->mod_pairs[k].i_2=j;
			cranck_data->mod_pairs[k].dist=-1.0;
			k++;
		}
	}
	cranck_data->N_pairs=k;
}

/**
 * Moves the ends of the protein backbone. How many residues from both ends
 * are affected is determined by #max_move_size.
 * @note The maximum angle is 0.1 radians (hardcoded)
 *
 * @param[in,out]   *pivot_data       Move data
 * @param[in,out]   *p                Modified protein
 * @param[in]       *rng_r            GSL random number generator state
 * @param[in]        max_move_size    Number of residues from each end to be affected by move
 *
 * @return \c void
 */
void CATMV_pivot	(mc_move_data *pivot_data, cat_prot *p, gsl_rng *rng_r, int max_move_size)
{
	int len,k;
	double v_1[3];
	double angle;
	int N_res=p->n_res;
	int verse;
	int type;
	int c,start;

    if (2*max_move_size > p->n_res)
    {
        error("Pivot move can not move more residues then number of residues in protein!", __FILE__, __LINE__);
    }

    c=gsl_rng_uniform_int (rng_r, 2*(max_move_size)); // select residue index in interval [0,(2*max_move_size)-1)
    if (c > max_move_size-1)
    {
        c=p->n_res-1-(c-max_move_size); // -1 since p->n_res is number of residues not number of maximal residue index
        verse=1;
    } else
    {
        verse=-1;
    }
	type=(gsl_rng_uniform(rng_r)>0.5);
	angle= 0.2*M_PI*(-0.5+gsl_rng_uniform(rng_r));
	if (verse>0) //rotate following atoms
	{
		if (type==0) //N-CA
		{
			//p->phi[c]+=angle;
			p->phi[c]= gsl_sf_angle_restrict_symm(p->phi[c]+angle);
			for(k=0;k<3;k++) { v_1[k] = p->CA[c][k] -	p->N [c][k]; }
			start=c*p->n_atom_per_res+ATOM_C;
			len=p->n_atoms-start;
			rotate(&p->coord[start],len,p->CA[c],v_1,angle);
		} else {
			//start from CA-C
			//p->psi[c]+=angle;
			p->psi[c]= gsl_sf_angle_restrict_symm(p->psi[c]+angle);
			for(k=0;k<3;k++) { v_1[k] = p->C[c][k] -	p->CA [c][k]; }
			start=c*p->n_atom_per_res+ATOM_O;
			len=p->n_atoms-start;
			if(p->CB!=NULL) {
				double CB[3]={p->CB[c][0],p->CB[c][1],p->CB[c][2]};
				rotate(&p->coord[start],len,p->C[c],v_1,angle);
				p->CB[c][0]=CB[0]; p->CB[c][1]=CB[1]; p->CB[c][2]=CB[2];
			} else {
				rotate(&p->coord[start],len,p->C[c],v_1,angle);
			}
		}
	} else {
		//rotate residues with index <c
		if (type==0) {
			//start from N-CA
			//p->phi[c]+=angle;
			p->phi[c]= gsl_sf_angle_restrict_symm(p->phi[c]+angle);
			for(k=0;k<3;k++) { v_1[k] = p->N[c][k] -	p->CA [c][k]; }
			len=c*p->n_atom_per_res+ATOM_N;//N excluded
			rotate(p->coord,len,p->N[c],v_1,angle);
		} else {
			//start from CA-C
			//p->psi[c]+=angle;
			p->psi[c]= gsl_sf_angle_restrict_symm(p->psi[c]+angle);
			for(k=0;k<3;k++) { v_1[k] = p->CA[c][k] -	p->C [c][k]; }
			len=c*p->n_atom_per_res+ATOM_CA; //CA excluded
			if(p->CB!=NULL) {
				rotate(p->coord,len,p->CA[c],v_1,angle);
				rotate(&p->CB[c],1,p->CA[c],v_1,angle);
			} else {
				rotate(p->coord,len,p->CA[c],v_1,angle);
			}
		}
	}
	//Store data on modified residues and residue pairs
	if (verse > 0) {
		pivot_data->N_moved=N_res-c;
		k=0;
		for(int i=0;i<pivot_data->N_moved;i++) {
			pivot_data->moved_res[i]=c+i;
			for(int j=0;j<=N_res-pivot_data->N_moved;j++) {
				pivot_data->mod_pairs[k].i_1=c+i;
				pivot_data->mod_pairs[k].i_2=j;
				pivot_data->mod_pairs[k].dist=-1.0;
				k++;
			}
		}
		pivot_data->N_pairs=k;
	} else {
		pivot_data->N_moved=c+1;
		k=0;
		for(int i=0;i<pivot_data->N_moved;i++) {
			//pivot_data->moved_res[i]=c+verse*i;
			pivot_data->moved_res[i]=i;
			for(int j=0;j<=N_res-pivot_data->N_moved;j++) {
				pivot_data->mod_pairs[k].i_1=i;
				pivot_data->mod_pairs[k].i_2=c+j;
				pivot_data->mod_pairs[k].dist=-1.0;
				k++;
			}
		}
		pivot_data->N_pairs=k;
	}
}

/**
 *
 * @param[in,out]   *ra_data          Move data
 * @param[in,out]   *p                Modified protein
 * @param[in]       *rng_r            GSL random number generator state
 * @param[in]        sigma            Sigma of normal distribution from which actual move length along tangent space is selected
 *
 * @return \c void
 */
void CATMV_concerted_rot(mc_move_data *ra_data, cat_prot *p, gsl_rng * rng_r, double sigma)
{
	// ra_data function only modifie which residues have moved so which pairs should be recalculated for E-change

	int
		k,                           // used how many and which pairs have to be modified
		start,                       // specife residue index from which concentrated rotation is performed
		error;                       // variable for handling error codes

	double
		w_ip1;                       // variable used in constraining dihedral angles in the range (-\pi,\pi]

	cr_input_data
		bb_in,                       // values used to store backbone initial configuration of 3 consecutive residues in format for concerted rotation
		bb_out;                      // modified values for backbone 3 consecutive residius generated by concerted rotation


	start= 1 + gsl_rng_uniform_int (rng_r, p->n_res - 4 );          // select random residue from which we start concerted rotation move

	bb_in = set_cr_params_phi(p,start);                             // allocate and initialize backbone parameters for 3 consecutive residues starting at residue with index "start"

	alloc_cr_input_data(&bb_out);                                   // allocate

	memcpy_cr_input_data(&bb_out, &bb_in);

	//----------CONCERTED ROTATION------------
	error = random_rot (bb_out, bb_in,rng_r, sigma);
	if(error!=GSL_SUCCESS) {
		free_cr_input_data(&bb_in);
		free_cr_input_data(&bb_out);
		return;
	}
	double w_ip2;
	for(int i = 0; i < 3; i++) {
		w_ip1 = gsl_vector_get(bb_out.dihed_angles, i*2)-M_PI;
		w_ip2 = gsl_vector_get(bb_out.dihed_angles, i*2+1);
		gsl_vector_set(bb_out.dihed_angles, 2*i, gsl_sf_angle_restrict_symm(w_ip1));      //These routines force the angle theta to lie in the range (-\pi,\pi]
		gsl_vector_set(bb_out.dihed_angles, 2*i+1, gsl_sf_angle_restrict_symm(w_ip2));      //These routines force the angle theta to lie in the range (-\pi,\pi]
	}
	// Rebuild the peptide chain base on set of dihedrals from concerted rotations
	for (int i = 0; i<3; i++) {
		CAT_add_peptide(
				p,
				start+i,
				gsl_vector_get(bb_out.dihed_angles,i*2),
				gsl_vector_get(bb_out.bend_angles, i*2),
				gsl_vector_get(bb_out.dihed_angles,(i*2)+1)
				);
	}
	free_cr_input_data(&bb_in);
	free_cr_input_data(&bb_out);
	// recalculate dihedral angles for whole protein
	for(int i = start-1; i < start+4; i++){
		if(i>0) {
			p->phi[i]=dihedralangle_ABCD(p->C[i-1],p->N[i],p->CA[i],p->C[i]);
		}
		if(i<p->n_res-1)
		{
			p->psi[i]=dihedralangle_ABCD(p->N[i],p->CA[i],p->C[i],p->N[i+1]);
		}
	}
	// Recalcuate the  positions of CB atoms if used
	if(p->CB!=NULL) {
		for(int i=start;i<start+4;i++) {
			CAT_insert_cbeta(p,i,CAT_Cb_AZIMUTH,CAT_Rbond_CCb);
		}
	}

	//save the data about moved pairs
	k=0;
	ra_data->N_moved=4;
	for(int i=0;i<ra_data->N_moved;i++) {
		ra_data->moved_res[i]=start+i;
		for(int j=0;j<start;j++) {
			ra_data->mod_pairs[k].i_1=start+i;
			ra_data->mod_pairs[k].i_2=j;
			ra_data->mod_pairs[k].dist=-1.0;
			k++;
		}
		for(int j=start+i;j<p->n_res;j++) {
			ra_data->mod_pairs[k].i_1=start+i;
			ra_data->mod_pairs[k].i_2=j;
			ra_data->mod_pairs[k].dist=-1.0;
			k++;
		}
	}
	ra_data->N_pairs=k;

	return;
}

/**
 *
 * @param[in] 	bb_in  starting configuration (following D-H convention) 
 * @param[in]   angle  free variable
 * @param[in]   delta_move size of the proposed displacement along the tangent space in bb_in (in radians).
 * @param[out]  bb_out new configuration (following D-H convention)
 *
 * @return GSL error code with the status of the computation. Returns GSL_SUCCESS
 * if a valid bb_out was found, some error code otherwise.
 * @note A rescaling of the step depending on the variable type (length or angle) 
 * should be considered in general, as specified in <a href="http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0118342">DOI:10.1371/journal.pone.0118342</a>. Here we do not, as we only move angles.
 */
int solver(cr_input_data bb_out, cr_input_data bb_in, int angle,double delta_mov)
{
	double T_ar[7];
	gsl_vector_view Tv=gsl_vector_view_array(T_ar,7);
	gsl_vector * T = &Tv.vector;
	gsl_matrix * B;
	//solver
  int status;
  size_t  iter = 0;
  const size_t n = 6; // dimension of multiroot solver ... how many equations are there
	const gsl_multiroot_fsolver_type *Ts;
  gsl_multiroot_fsolver *s;
	double x_ar[6];
	gsl_vector_view x_v=gsl_vector_view_array(x_ar,6);
  gsl_vector *x = &x_v.vector;
	//function parameters
  struct rparams p;
	//choose an angle, get tangential vector to manifold
	switch (angle) {
		case 0:
			status=TmT7_t1(T,bb_in);
			break;
		case 1:
			status=TmT7_t2(T,bb_in);
			break;
		case 2:
			status=TmT7_t3(T,bb_in);
			break;
		case 3:
			status=TmT7_t4(T,bb_in);
			break;
		case 4:
			status=TmT7_t5(T,bb_in);
			break;
		case 5:
			status=TmT7_t6(T,bb_in);
			break;
		case 6:
			status=TmT7_t7(T,bb_in);
			break;
	}
	if(status!=0) {
		//printf("Jacobian not invertible. dice: %d Error code %d\n",angle,status);
		return GSL_FAILURE; //again, I know it is not gsl's fault..
	}
  //printf ("angle %d status = %s\n",angle, gsl_strerror (status));
	B=gram_schmidt (T);
	//parameters
	//printf("AJJJJJ delta_mov=%lf\n",delta_mov);
	p.s=delta_mov;
	p.B=B;
	//p.wspace=gsl_vector_alloc(7);
	double ws_ar[7];
	gsl_vector_view ws_v=gsl_vector_view_array(ws_ar,7);
	p.wspace=&ws_v.vector;
	//
	p.bb_out=bb_out;
	p.bb_in	=bb_in;
	//
  gsl_multiroot_function f = {&T7_solver, n, &p};
	gsl_vector_set_zero(x);
	//solver inizialization
  Ts = gsl_multiroot_fsolver_hybrid;
  s = gsl_multiroot_fsolver_alloc (Ts, 6);
  gsl_multiroot_fsolver_set (s, &f, x);

  //print_state (iter, s);
  do {
      iter++;
      status = gsl_multiroot_fsolver_iterate (s);
      //print_state (iter, s);
			/* check if solver is stuck */
      if (status) {
				//printf ("LEAVING---status = %s. iter %d\n", gsl_strerror (status),iter);
        break;
			}
      status = gsl_multiroot_test_residual (s->f, 1e-12);
	} while (status == GSL_CONTINUE && iter < 100);
	//if(iter>=100) { printf ("---status = %s. angle %d\n", gsl_strerror (status),angle+1);}
/*
	//check parity preservation
	if(status==GSL_SUCCESS)
	{
		double T7_33_old,T7_33_new;
		struct expl_data old,new;
		for(int i=0;i<7;i++)
		{
			old.r[i]=gsl_vector_get(bb_in.r,i);
			old.c[i]=cos(gsl_vector_get(bb_in.dihed_angles,i));
			old.s[i]=sin(gsl_vector_get(bb_in.dihed_angles,i));
			old.ca[i]=cos(gsl_vector_get(bb_in.bend_angles,i));
			old.sa[i]=cos(gsl_vector_get(bb_in.bend_angles,i));
			new.r[i]=gsl_vector_get(bb_out.r,i);
			new.c[i]=cos(gsl_vector_get(bb_out.dihed_angles, i));
			new.s[i]=sin(gsl_vector_get(bb_out.dihed_angles, i));
			new.ca[i]=cos(gsl_vector_get(bb_out.bend_angles,i));
			new.sa[i]=cos(gsl_vector_get(bb_out.bend_angles,i));
		}
		T7_33_old=get_T7_33(old);
		T7_33_new=get_T7_33(new);
		if(GSL_SIGN(T7_33_old) != GSL_SIGN(T7_33_new))
		{
			//gsl_vector_memcpy(xi_new,xi_old);
			status = GSL_FAILURE; //YEP I know it is not gsl fault.
			fprintf(stderr,"Parity violated!!\n");
		} //devo controllare che non overlappi roba della gsl...
		//if parity is violated f and v are now different, by a sign.
	}
*/
  //gsl_vector_free (p.wspace);
  gsl_multiroot_fsolver_free (s);
  gsl_matrix_free (B);
	return status;
}

/**
 *
 * @param[in,out]   *x                Solution vector
 * @param[in,out]   *p                Parameters
 * @param[in]       *f                Right side vector
 *
 * @return GLS error value
 */
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
	for(i=0;i<6;i++)
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
	gsl_vector_view  w_ = gsl_vector_subvector (w, 1, 6);
	T7( f, bb_out );
	T7( &w_.vector, bb_in);
	gsl_vector_sub(f,&w_.vector);
	return error;
}

/**
 * @param[in]        iter             Iteration number
 * @param[in]       *s                GLS multi root solver
 *
 * @return \c void
 */
void print_state (int iter, gsl_multiroot_fsolver * s)
{
  printf ("iter = %3u x = % .3f % .3f % .3f % .3f % .3f % .3f"
          "f(x) = % .3f % .3f % .3f % .3f% .3f % .3f\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->x, 2),
          gsl_vector_get (s->x, 3),
          gsl_vector_get (s->x, 4),
          gsl_vector_get (s->x, 5),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1),
          gsl_vector_get (s->f, 2),
          gsl_vector_get (s->f, 3),
          gsl_vector_get (s->f, 4),
          gsl_vector_get (s->f, 5));
}

/**
 * @param[in]        a                Index to which angle determinant is calculated
 * @param[in]        bb_in            Current configuration of protein backbone in concerted rotation format
 *
 * @return determinant
 */
double jac_Det(int a, cr_input_data bb_in)
{
	int error, sign;
	double J;
	//int i;
	//double x,J;
	//double c[7],s[7];
	gsl_matrix *M;
	gsl_permutation * p = gsl_permutation_alloc (6);
	switch (a) {
		case 0:
			M=jac_t1(bb_in);
			break;
		case 1:
			M=jac_t2(bb_in);
			break;
		case 2:
			M=jac_t3(bb_in);
			break;
		case 3:
			M=jac_t4(bb_in);
			break;
		case 4:
			M=jac_t5(bb_in);
			break;
		case 5:
			M=jac_t6(bb_in);
			break;
		case 6:
			M=jac_t7(bb_in);
			break;
	}
	error=gsl_linalg_LU_decomp(M,p,&sign);
	J=gsl_linalg_LU_det(M,sign);
	gsl_matrix_free(M);
	gsl_permutation_free(p);
	return J;
}

/**
 * @note Parameters set starts with phi angle of first residue defined by c value
 *
 * @param[in]       *p                Backbone reprezentation in #cat_prot format
 * @param[in]        c                Index of first residue to generate input parameters for concerted rotation
 *
 * @return parameters for concerted rotation move
 */
cr_input_data set_cr_params_phi( cat_prot *p, int c)
{
	int
		i,
		j;

	double
		r,
		z_i[3],
		z_im1[3],
		s[3],
		d[3];

		cr_input_data
			cr_in;

		alloc_cr_input_data(&cr_in);
		//params - PHI matrix
		for(i=0;i<4;i++)
		{
			j=i*2;
			//r
			z_i[0]=p->CA[c+i][0]-p->N[c+i][0];
			z_i[1]=p->CA[c+i][1]-p->N[c+i][1];
			z_i[2]=p->CA[c+i][2]-p->N[c+i][2];
			//Questo  e` z_ip1...ma tant'e`...
			z_im1[0]=p->C[c+i][0]-p->CA[c+i][0];
			z_im1[1]=p->C[c+i][1]-p->CA[c+i][1];
			z_im1[2]=p->C[c+i][2]-p->CA[c+i][2];
			gsl_vector_set(cr_in.r,j,norm_d(&z_i[0],3));
			//bending angles
			normalize_d(&z_i[0],3);
			normalize_d(&z_im1[0],3);
			gsl_vector_set(cr_in.bend_angles,j,acos(scal_d(z_i,z_im1,3)));
			gsl_vector_set(cr_in.d,j,0);
			//gsl_vector_set(cr_in.bend_angles,j,acos(scal_d(z_i,z_im1,3)));
			//gsl_vector_set(cr_in.bend_angles,j,M_PI-1.20428);
			//dihedrals
			gsl_vector_set(cr_in.dihed_angles,j,p->phi[c+i]+M_PI);
		}
		//params - PSI matrix
		for(i=0;i<3;i++)
		{
			j=2*i+1;
			//r
			z_i[0]=p->C[c+i][0]-p->CA[c+i][0];
			z_i[1]=p->C[c+i][1]-p->CA[c+i][1];
			z_i[2]=p->C[c+i][2]-p->CA[c+i][2];
			//again, z_ip1..
			z_im1[0]=p->CA[c+i+1][0]-p->N[c+i+1][0];
			z_im1[1]=p->CA[c+i+1][1]-p->N[c+i+1][1];
			z_im1[2]=p->CA[c+i+1][2]-p->N[c+i+1][2];
			s[0]=p->N[c+i+1][0]-p->CA[c+i][0];
			s[1]=p->N[c+i+1][1]-p->CA[c+i][1];
			s[2]=p->N[c+i+1][2]-p->CA[c+i][2];
			normalize_d(&z_i[0],3);
			normalize_d(&z_im1[0],3);
			r=scal_d(&s[0],&z_i[0],3); // project s vector to Ca-C vector
			gsl_vector_set(cr_in.r,j,r);
			//d
			d[0]=s[0]-r*z_i[0]; //distance from end of S projected on z_i to N[c+i+1]
			d[1]=s[1]-r*z_i[1];
			d[2]=s[2]-r*z_i[2];
			gsl_vector_set(cr_in.d,j,norm_d(&d[0],3));
			//normalize_d(&s[0],3);
			//bending angles
			gsl_vector_set(cr_in.bend_angles,j,acos(scal_d(z_i,z_im1,3)));
			//gsl_vector_set(cr_in.bend_angles,j,acos(scal_d(z_i,z_im1,3)));
			//gsl_vector_set(cr_in.bend_angles,j,0.10995574287);
			//dihedrals
			gsl_vector_set(cr_in.dihed_angles,j,p->psi[c+i]);
		}
		return cr_in;
}

/**
 * @note Parameters set starts with psi angle of first residue defined by c value
 *
 * @param[in]       *p                Backbone reprezentation in #cat_prot format
 * @param[in]        c                Index of first residue to generate input parameters for concerted rotation
 *
 * @return parameters for concerted rotation move
 */
cr_input_data set_cr_params_psi( cat_prot *p, int c)
{
	int
		i,
		j;

	double
		r,
		z_i[3],
		z_im1[3],
		s[3],
		d[3];

		cr_input_data
			cr_in;

		alloc_cr_input_data(&cr_in);
		//params - PHI matrix
		for(i=0;i<3;i++)
		{
			j=i*2+1;
			//r
			z_i[0]=p->CA[c+i][0]-p->N[c+i][0];
			z_i[1]=p->CA[c+i][1]-p->N[c+i][1];
			z_i[2]=p->CA[c+i][2]-p->N[c+i][2];
			//Questo  e` z_ip1...ma tant'e`...
			z_im1[0]=p->C[c+i][0]-p->CA[c+i][0];
			z_im1[1]=p->C[c+i][1]-p->CA[c+i][1];
			z_im1[2]=p->C[c+i][2]-p->CA[c+i][2];
			gsl_vector_set(cr_in.r,j,norm_d(&z_i[0],3));
			//bending angles
			normalize_d(&z_i[0],3);
			normalize_d(&z_im1[0],3);
			gsl_vector_set(cr_in.bend_angles,j,acos(scal_d(z_i,z_im1,3)));
			//gsl_vector_set(cr_in.bend_angles,j,acos(scal_d(z_i,z_im1,3)));
			//gsl_vector_set(cr_in.bend_angles,j,M_PI-1.20428);
			//dihedrals
			gsl_vector_set(cr_in.d,j,0);
			gsl_vector_set(cr_in.dihed_angles,j,p->phi[c+i]+M_PI);
		}
		//params - PSI matrix
		for(i=0;i<4;i++)
		{
			j=2*i;
			//r
			z_i[0]=p->C[c+i][0]-p->CA[c+i][0];
			z_i[1]=p->C[c+i][1]-p->CA[c+i][1];
			z_i[2]=p->C[c+i][2]-p->CA[c+i][2];
			//again, z_ip1..
			z_im1[0]=p->CA[c+i+1][0]-p->N[c+i+1][0];
			z_im1[1]=p->CA[c+i+1][1]-p->N[c+i+1][1];
			z_im1[2]=p->CA[c+i+1][2]-p->N[c+i+1][2];
			s[0]=p->N[c+i+1][0]-p->CA[c+i][0];
			s[1]=p->N[c+i+1][1]-p->CA[c+i][1];
			s[2]=p->N[c+i+1][2]-p->CA[c+i][2];
			normalize_d(&z_i[0],3);
			normalize_d(&z_im1[0],3);
			r=scal_d(&s[0],&z_i[0],3); // project s vector to Ca-C vector
			gsl_vector_set(cr_in.r,j,r);
			//d
			d[0]=s[0]-r*z_i[0]; //distance from end of S projected on z_i to N[c+i+1]
			d[1]=s[1]-r*z_i[1];
			d[2]=s[2]-r*z_i[2];
			gsl_vector_set(cr_in.d,j,norm_d(&d[0],3));
			//normalize_d(&s[0],3);
			//bending angles
			gsl_vector_set(cr_in.bend_angles,j,acos(scal_d(z_i,z_im1,3)));
			//gsl_vector_set(cr_in.bend_angles,j,acos(scal_d(z_i,z_im1,3)));
			//gsl_vector_set(cr_in.bend_angles,j,0.10995574287);
			//dihedrals
			gsl_vector_set(cr_in.dihed_angles,j,p->psi[c+i]);
		}
		return cr_in;
}

/**
 * Takes the initial backbone configuration generated either via #set_cr_params_phi or #set_cr_params_psi and width of normal distribution sigma from where
 * move size is selected and return modified backbone configuration in concerted rotation format.
 *
 * @param[in,out]    bb_out           new D-H backbone configuration 
 * @param[in,out]    bb_in            starting D-H backbone configuration 
 * @param[in,out]   *rng_r            State of GLS random generator
 * @param[in]        sigma            Width of normal distribution used to select the move length along the tangent space
 *
 * @return GLS error code
 */
int random_rot(cr_input_data bb_out, cr_input_data bb_in, gsl_rng *rng_r, double sigma)
{
	int error,a,i,m;
	int angles_idx[7]={0,1,2,3,4,5,6};
	double ds;
	//Try all 7 rotation angles.
	m=7;
	ds=gsl_ran_gaussian(rng_r,sigma);
	error=1;
	do{
		i=gsl_rng_uniform_int (rng_r, m);
		a=angles_idx[i];
		error=solver(bb_out,bb_in,a,ds);
		if(error) {
			m--;
			angles_idx[i]=angles_idx[m];
		}
	}while(m && error!=GSL_SUCCESS);
	//Work out the correct acceptance
	double ds1;
	double xi_new_ar[7];
	double delta_xi_ar[7];
	gsl_vector_view T_xi_new_v =gsl_vector_view_array(xi_new_ar,7);
	gsl_vector *T_xi_new=&T_xi_new_v.vector;
	gsl_vector_view delta_xi_v =gsl_vector_view_array(delta_xi_ar,7);
	gsl_vector *delta_xi=&delta_xi_v.vector;
	//To be improved..
	if(error==GSL_SUCCESS) {
		double Jo,Jn;
		Jo=fabs(jac_Det(a,bb_in));
		m=7;
		int status=1;
		int angles_idx_B[7]={0,1,2,3,4,5,6};
		do{
			i=gsl_rng_uniform_int (rng_r, m);
			a=angles_idx_B[i];
			switch (a) {
				case 0:
					status=TmT7_t1(T_xi_new,bb_out);
					break;
				case 1:
					status=TmT7_t2(T_xi_new,bb_out);
					break;
				case 2:
					status=TmT7_t3(T_xi_new,bb_out);
					break;
				case 3:
					status=TmT7_t4(T_xi_new,bb_out);
					break;
				case 4:
					status=TmT7_t5(T_xi_new,bb_out);
					break;
				case 5:
					status=TmT7_t6(T_xi_new,bb_out);
					break;
				case 6:
					status=TmT7_t7(T_xi_new,bb_out);
					break;
			}
			if(status) {
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
		if(error==GSL_SUCCESS) {
			double w1,w2;
			for(i=0;i<7;i++) {
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

