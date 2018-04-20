
#include "geom_prop.h"


/**
 * @brief Function generate inertia tensor from point coordinates and their masses
 *
 *
 * @param[in]     **coord            Set of points.
 * @param[in]      *mass             Set of masses for each point in **coord\.
 * @param[in]       N                Number of points/masses in **coord or *mass arrays.
 * @param[out]     *inrt             GSL matrix containing inertia tensor.
 *
 * @return \c void
 */
void inertia_tensor (  double **coord, double *mass,int N, gsl_matrix* inrt )
{
  int i,j;
  double cm[3];
	double v[3];
  double **cm_coord;
  double xx, yy, zz;
  double xy, xz, yz;
  double M=0; // Mass of whole object
  xx=yy=zz=0.0;
  xy=xz=yz=0.0;
	cm[0]=cm[1]=cm[2]=0;
	for(i=0;i<N;i++)
	{
		for(j=0;j<3;j++) { cm[j]+=mass[i]*coord[i][j]; }
		M+=mass[i];
	}
	for(j=0;j<3;j++){ cm[j]/= N; }
  for(i=0;i<N;i++)
	{
		for(j=0;j<3;j++){v[j]=coord[i][j]-cm[j];}
		xx+= mass[i]*(v[1]*v[1]+v[2]*v[2]);
		yy+= mass[i]*(v[0]*v[0]+v[2]*v[2]);
		zz+= mass[i]*(v[0]*v[0]+v[1]*v[1]);
		xy-= mass[i]*v[0]*v[1];
		xz-= mass[i]*v[0]*v[2];
		yz-= mass[i]*v[1]*v[2];
  }
  gsl_matrix_set(inrt,0,0,(xx/(2*M*M)));
  gsl_matrix_set(inrt,1,1,(yy/(2*M*M)));
  gsl_matrix_set(inrt,2,2,(zz/(2*M*M)));
  gsl_matrix_set(inrt,0,1,(xy/(2*M*M)));
  gsl_matrix_set(inrt,1,0,(xy/(2*M*M)));
  gsl_matrix_set(inrt,0,2,(xz/(2*M*M)));
  gsl_matrix_set(inrt,2,0,(xz/(2*M*M)));
  gsl_matrix_set(inrt,1,2,(yz/(2*M*M)));
  gsl_matrix_set(inrt,2,1,(yz/(2*M*M)));
}

/**
 * @brief Function generate #SYMM_TENS tensor structure from point coordinates and their masses
 *
 * Generation of inertia tensor inside function is done by inertia_tensor()\.
 * Eigenvectors and values are then calculated from inertia tensor.
 *
 * @param[in]     **coord            Set of points.
 * @param[in]      *mass             Set of masses for each point in **coord.
 * @param[in]       N                Number of points/masses in **coord or *mass arrays.
 *
 * @return #SYMM_TENS structure
 */
SYMM_TENS * inertia_axes( double **coord, double *mass, int N)
{
	SYMM_TENS *inertia;
	if((inertia=(SYMM_TENS*)malloc(sizeof(SYMM_TENS)))==NULL)
	{
		failed("failed allocation in inertia_axes");
	}
	inertia->size=3;
	inertia->matrix	= gsl_matrix_alloc(3,3);
	inertia->evecs 	= gsl_matrix_alloc(3,3);
	inertia->evals 	= gsl_vector_alloc(3);
	gsl_eigen_symmv_workspace *sym_w_spc = gsl_eigen_symmv_alloc(3);
	inertia_tensor	(coord,mass,N,inertia->matrix);
	gsl_eigen_symmv (inertia->matrix, inertia->evals,inertia->evecs, sym_w_spc);
	gsl_eigen_symmv_sort (inertia->evals, inertia->evecs, GSL_EIGEN_SORT_VAL_DESC);
	gsl_eigen_symmv_free(sym_w_spc);
	return inertia;
}

/**
 * @brief Function calculate asphericity
 *
 * @todo Add more info.
 *
 * @param[in]      *semiaxis        3D vector.
 *
 * @return asphericity coeficient
 */
double asphericity(gsl_vector *semiaxis){
  double asph;
  double a,b,c;
  gsl_sort_vector(semiaxis);
  a=gsl_vector_get(semiaxis,2);
  b=gsl_vector_get(semiaxis,1);
  c=gsl_vector_get(semiaxis,0);

  asph=0.5*((a-b)*(a-b)+(a-c)*(a-c)+(b-c)*(b-c))/pow((a+b+c),2);
  return asph;
}

/**
 * @brief Function calculate prolateness
 *
 * @todo Add more info.
 *
 * @param[in]      *semiaxis        3D vector.
 *
 * @return prolateness coeficient
 */
double prolateness(gsl_vector *semiaxis){
  double prol;
  double den;
  double a,b,c;
  gsl_sort_vector(semiaxis);
  a=gsl_vector_get(semiaxis,2);
  b=gsl_vector_get(semiaxis,1);
  c=gsl_vector_get(semiaxis,0);
  den=(a*a+b*b+c*c-a*b-a*c-b*c);
  den=pow(den,3./2.);
  prol=0.5*(2*a-b-c)*(2*b-a-c)*(2*c-a-b)/den;
  return prol;
}

/**
 * @brief Function calculate radius of gyration
 *
 *
 * @param[in]     **coord        N 3D points.
 * @param[in]       N            Number of points.
 *
 * @return radius of gyration
 */
double gyration_radius(double **coord,int N)
{
  int i;
  double cm[3], r2;
  r2=0;
  center_of_mass(N,coord,cm);
  for(i=0;i<N;i++)
    r2+=pow(dist_d(coord[i],cm,3),2);
  r2=r2/N;
  return sqrt(r2);
}

/** @brief Compute envelope tensor???
 */
/*
void XXX_tensor(double **coord, int N,  double *mass, gsl_matrix* inrt)
{
  int i,j;
  double cm[3];
  double **cm_coord;
  double xx, yy, zz;
  double xy, xz, yz;
  double M=0;
  xx=yy=zz=0.0;
  xy=xz=yz=0.0;
  center_of_mass(N,coord,cm);
  for(i=0;i<N;i++)
  {
    M+=mass[i];
  }
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      xx+=(coord[i][0]-coord[j][0])*(coord[i][0]-coord[j][0]);
      yy+=(coord[i][1]-coord[j][1])*(coord[i][1]-coord[j][1]);
      zz+=(coord[i][2]-coord[j][2])*(coord[i][2]-coord[j][2]);
      xy+=(coord[i][0]-coord[j][0])*(coord[i][1]-coord[j][1]);
      xz+=(coord[i][0]-coord[j][0])*(coord[i][2]-coord[j][2]);
      yz+=(coord[i][1]-coord[j][1])*(coord[i][2]-coord[j][2]);
    }
  }
  gsl_matrix_set(inrt,0,0,(xx/(2*M*M)));
  gsl_matrix_set(inrt,1,1,(yy/(2*M*M)));
  gsl_matrix_set(inrt,2,2,(zz/(2*M*M)));
  gsl_matrix_set(inrt,0,1,(xy/(2*M*M)));
  gsl_matrix_set(inrt,1,0,(xy/(2*M*M)));
  gsl_matrix_set(inrt,0,2,(xz/(2*M*M)));
  gsl_matrix_set(inrt,2,0,(xz/(2*M*M)));
  gsl_matrix_set(inrt,1,2,(yz/(2*M*M)));
  gsl_matrix_set(inrt,2,1,(yz/(2*M*M)));
  
}
*/


/**
 * @brief Function rotate set of N points along given center and rotation axis by given angle
 *
 *
 * @param[in,out] **coord        N 3D points to be rotated.
 * @param[in]       N            Number of points.
 * @param[in]       center       3D point used as center of rotation.
 * @param[in]       axis         Vector around which rotation is performed.
 * @param[in]       theta        Angle by which set of points is rotated.
 *
 * @return \c void
 */
void rotate ( double **coord, int N,double center[3],double axis[3], double theta)
{
  int i;
  double CsTh, SnTh,OnemCsTh;
	double nx,ny,nz;
  //double N[3];
	if (N == 0) {return;}
  gsl_vector *v = gsl_vector_alloc (3);
  gsl_vector *w = gsl_vector_alloc (3);
  gsl_matrix *R = gsl_matrix_alloc(3,3);

  if ( coord == NULL || coord[0]==NULL)
  {
    failed("tryng to rotate an empty protein pointer!\n");
  }

  normalize_d(axis,3);
  nx = axis[0]; ny = axis[1]; nz = axis[2]; // I know it is dumb, just reusing code..
  CsTh = cos(theta);
  SnTh = sin(theta);
  OnemCsTh = 1 - CsTh;
  //random versor for PIVOT rotation |n| = 1
  //rotation matrix, from quaternions multiplication
  //first row
  gsl_matrix_set ( R,0,0,CsTh + nx*nx* OnemCsTh);
  gsl_matrix_set ( R,0,1,nx*ny*OnemCsTh - nz*SnTh);
  gsl_matrix_set ( R,0,2,nx*nz*OnemCsTh + ny*SnTh);
  //second row
  gsl_matrix_set ( R,1,0,nx*ny*OnemCsTh +nz*SnTh);
  gsl_matrix_set ( R,1,1,CsTh + ny*ny*OnemCsTh);
  gsl_matrix_set ( R,1,2,ny*nz*OnemCsTh - nx*SnTh);
  //third row
  gsl_matrix_set ( R,2,0,nx*nz*OnemCsTh - ny*SnTh);
  gsl_matrix_set ( R,2,1,ny*nz*OnemCsTh + nx*SnTh);
  gsl_matrix_set ( R,2,2,CsTh + nz*nz*OnemCsTh);
  for ( i = N-1; i >=0 ; i-- )
  {
    coord[i][0] -=center[0];
    coord[i][1] -=center[1];
    coord[i][2] -=center[2];
  }
  for ( i = N-1; i >=0 ; i-- )
  {
    gsl_vector_set ( v, 0, coord[i][0]);
    gsl_vector_set ( v, 1, coord[i][1]);
    gsl_vector_set ( v, 2, coord[i][2]);
    gsl_blas_dgemv ( CblasNoTrans, 1, R, v, 0,w);

    coord[i][0] = gsl_vector_get(w,0)+center[0];
    coord[i][1] = gsl_vector_get(w,1)+center[1];
    coord[i][2] = gsl_vector_get(w,2)+center[2];
  }
	gsl_vector_free(v);
	gsl_vector_free(w);
	gsl_matrix_free(R);
}

/**
 * @brief Function create orthonormal basis with respect to 3D vector
 *
 * Metod use Gram-Schmidt orthogonalization.
 *
 * @todo Rename fnction to gram_schmidt ... might be an typo.
 *
 * @param[in]      *v            3D vector from which base is calculated.
 *
 * @return Orhonormal basis in 3D with *v as one of coordinates
 */
gsl_matrix * graham_schmidt ( gsl_vector *v)
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
