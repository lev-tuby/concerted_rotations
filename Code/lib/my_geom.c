/**
 * @file
 * @brief Source file for all functions for work vectors and geometric object. Adapted from the <a href="https://github.com/luca-tubiana/KymoKnot"> KymoKnot github repository </a> by Luca Tubiana, Guido Polles, Enzo Orlandini and Cristian Micheletti.
 * Vectors are implemented as arrays without boundary checks.
 */

#include "my_geom.h"
#include "messages.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_blas.h>

/** @brief Machine floating point precision. */
#ifndef	EPSIL
#define EPSIL  0.000000000001
#endif

/************ GEOMETRIC ROUTINES ************/
/*******************************/
/**
 *
 * @param[in]      *a                First 3D vector.
 * @param[in]      *b                Second 3D vector.
 * @param[out]     *c                3D vector \f$\vec{c}=\vec{a}\times \vec{b}\f$.
 *
 * @return \c void
 */
void vecprod_d (const double *a, const double *b, double *c)
{

  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}
/*******************************/

/**
 *
 * @note Assumes both vectors are of dimensin dim. Otherwise it might cause a SEGFAULT.
 *
 * @param[in]      *a                First dim D dimensional vector.
 * @param[in]      *b                Second dim D dimensional vector.
 * @param[in]       dim              Dimension of vectors.
 *
 * @return scalar product \f$\sum_i=1^{dim} a_i b_i\f$.
 */
double scal_d (const double *a, const double *b, int dim)
{

	int i;
	double temp;
  temp = 0.0;
	switch(dim) {
		case 1:
			temp=a[0]*b[0];
			break;
		case 2:
			temp=a[0]*b[0]+a[1]*b[1];
			break;
		case 3:
			temp=a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
			break;
		default:
			for (i = 0; i < dim; i++) {
				temp += a[i] * b[i];
			}
			break;
	}
  return (temp);
}
/*******************************/


/**
 *
 * @note Assumes both vectors are of dimensin dim. Otherwise it might cause a SEGFAULT.
 *
 * @param[in]      *a                First dim D dimensional vector.
 * @param[in]       dim              Dimension of vectors.
 *
 * @return length(norm) of vector *a
 */
double norm_d (const double *a, int dim)
{
  return (sqrt (scal_d (a, a, dim)));
}
/*******************************/

/*******************************/

/**
 *
 * @note Assumes both vectors are of dimensin dim. Otherwise it might cause a SEGFAULT.
 *
 * @param[in,out]  *a                dim D dimensional vector to be normalized.
 * @param[in]       dim              Dimension of vectors.
 *
 * @return length(norm) of vector *a
 */
void normalize_d (double *a, int dim)
{
	int i;
	double temp;

	temp = norm_d (a, dim);
	for (i = 0; i < dim; i++) {
		a[i] = a[i] / temp;
	}
}



/*******************************/
/**
 *
 * @note Assumes both vectors are of dimensin dim. Otherwise it might cause a SEGFAULT.
 *
 * @param[in]      *a                First dim D dimensional vector.
 * @param[in]      *b                Second dim D dimensional vector.
 * @param[in]       dim              Dimension of vectors.
 *
 * @return distance between points
 */
double dist_d (const double *a, const double *b, int dim)
{

	unsigned short i;
	double temp;

	temp = 0.0;
	for (i = 0; i < dim; i++) {
		temp += (a[i] - b[i]) * (a[i] - b[i]);
	}

	temp = sqrt (temp);
	return (temp);
}

/*******************************/
/**
 *
 * First segment is defined as vector p1 - p0 and second as p3 - p2.
 * @note All vectors have to be in 3D.
 *
 * @param[in]       p0               First segment beginning point.
 * @param[in]       p1               First segment ending point.
 * @param[in]       p2               Second segment beginning point.
 * @param[in]       p3               Second segment ending point.
 *
 * @return minimal distance two line segments
 */
double dist_segments(const double p0[3], const double p1[3], const double p2[3], const double p3[3])
{

  double d02[3],d32[3],d10[3], pa[3], pb[3], mua, mub, min_dist;
  int k;
  double d0232,d3210,d0210,d3232,d1010;

  for(k=0; k < 3; k++){
    d02[k]=p0[k]-p2[k];
    d32[k]=p3[k]-p2[k];
    d10[k]=p1[k]-p0[k];
  }

  d0232=0.0;
  d3210=0.0;
  d0210=0.0;
  d3232=0.0;
  d1010=0.0;

  for(k=0; k < 3; k++){
    d0232+=d02[k]*d32[k];
    d3210+=d32[k]*d10[k];
    d0210+=d02[k]*d10[k];
    d3232+=d32[k]*d32[k];
    d1010+=d10[k]*d10[k];
  }

  mua= (d0232*d3210-d0210*d3232)/(d1010*d3232-d3210*d3210);

  mub=(d0232+mua*d3210)/d3232;

  if (mua <0) mua=0;
  else if (mua >1) mua=1;

  if (mub <0) mub=0;
  else if (mub >1) mub=1;

  for(k=0; k < 3; k++){
    pa[k]=p0[k]+mua*(p1[k]-p0[k]);
    pb[k]=p2[k]+mub*(p3[k]-p2[k]);
  }

  min_dist=dist_d(pa,pb,3);

  return(min_dist);
}

/**
 *
 * Determines whether or not the line segment p1,p2
 * Intersects the 3 vertex facet bounded by pa,pb,pc
 * Return true/false and the intersection point p
 * The equation of the line is p = p1 + mu (p2 - p1)
 * The equation of the plane is a x + b y + c z + d = 0
 * n.x x + n.y y + n.z z + d = 0
 *
 * @param[in]       p1               Begining point of first segment.
 * @param[in]       p2               End point of first segment.
 * @param[in]       pa               First point defining vertex.
 * @param[in]       pb               Second point defining vertex.
 * @param[in]       pc               Third point defining vertex.
 * @param[out]     *p                If function return #TRUE point of intersection.
 *
 * @return #TRUE if segment intersect vertex otherwise #FALSE
 */
int LineFacet(const double *p1, const double *p2, const double *pa, const double *pb, const double *pc, double *p)
{
 /* new version 18/5/2009, CM and LT */
   double d;
   double bond[3];
   /*double a1,a2,a3;*/
   double denom,mu;
   double n[3];/*pa1[3],pa2[3],pa3[3];*/
   unsigned short k;
   double vp[3], v1[3], v2[3];
   double scal_vpv1, scal_v1v1, scal_v2v2, scal_v1v2, scal_vpv2;
   double x,y;
   double EPSIL2; /* tolerance for intersection checks */

   EPSIL2=0.0001; /* higher values make the criterion more stringent (but possibly discard
                     genuinely OK configurations. previously 0.0001*/

   for(k=0; k < 3; k++){
     v1[k]   = pa[k] - pb[k];
     v2[k]   = pc[k] - pb[k];
     bond[k] = p2[k] - p1[k];
   }

 /* Calculate the parameters for the plane */
   n[0] = v1[1]*v2[2] - v1[2]*v2[1];
   n[1] = v1[2]*v2[0] - v1[0]*v2[2];
   n[2] = v1[0]*v2[1] - v1[1]*v2[0];
   //n[0] = (pb[1] - pa[1])*(pc[2] - pa[2]) - (pb[2] - pa[2])*(pc[1] - pa[1]);
   //n[1] = (pb[2] - pa[2])*(pc[0] - pa[0]) - (pb[0] - pa[0])*(pc[2] - pa[2]);
   //n[2] = (pb[0] - pa[0])*(pc[1] - pa[1]) - (pb[1] - pa[1])*(pc[0] - pa[0]);
   normalize_d(n,3);
   d = - scal_d(n,pa,3);

/* Calculate the position on the line that intersects the plane */
   denom = n[0] * bond[0] + n[1] * bond[1] + n[2] * bond[2];
   if (fabs(denom) < EPSIL){         /* Line and plane don't intersect */
     return(FALSE);
   }

   mu = - (d + scal_d(n,p1,3))/ denom;


   if (mu <  (-EPSIL2) || mu > (1.0 + EPSIL2))   /* Intersection not along line segment */
     return(FALSE);

     /* Intersection point */
     for(k=0; k < 3; k++){
       p[k] = p1[k] + mu * bond[k];
       vp[k] = p[k]  - pb[k];
     }


     /* Determine whether or not the intersection point is bounded by pa,pb,pc */

   /* we write the location of the intersection point as a linear combination of two
      vectors delimiting the facet. The origin is taken as pb.
      The facet is the locus of points {a*v1+b*v2 | a,b >=0 & a+b <=1}. */

   scal_vpv1 = scal_d(vp,v1,3);
   scal_vpv2 = scal_d(vp,v2,3);
   scal_v1v2 = scal_d(v1,v2,3);
   scal_v1v1 = scal_d(v1,v1,3);
   scal_v2v2 = scal_d(v2,v2,3);

   denom =  scal_v1v1*scal_v2v2 - scal_v1v2*scal_v1v2;

   x= ((scal_v2v2)*(scal_vpv1) - (scal_v1v2)*(scal_vpv2))/denom;
   y= -((scal_v1v2)*(scal_vpv1) - (scal_v1v1)*(scal_vpv2))/denom;


   if ((x <  (-EPSIL2) ) || (y <  (-EPSIL2) ))   return(FALSE);
   if ( (x+y) > (1 + EPSIL2))   return(FALSE);

   return(TRUE);
}


/**
 *
 *
 * @param[in]       N                Number of of points (N<= number of coords!).
 * @param[in]     **coord            Array of 3D points.
 * @param[out]     *cm               Center of mass of N points in **coord.
 *
 * @return \c void
 */
void center_of_mass(int N, double ** coord, double *cm)
{
  int i;
  cm[0] = cm[1] = cm[2] = 0.0;
  for ( i = 0 ; i < N ; i++ ) {
    cm[0] += coord[i][0];
    cm[1] += coord[i][1];
    cm[2] += coord[i][2];
  }

  cm[0] /= N;
  cm[1] /= N;
  cm[2] /= N;
}

/**
 *
 *
 * @param[in]       N                Number of of points (N<= number of coords!).
 * @param[in]     **coord            Array of 3D points.
 * @param[out]      center           Center of mass of N points in **coord.
 *
 * @return distance of farthest point from CM
 * @todo speed up the computation by avoiding computing the sqrt N times.
 */
double find_radius ( int N ,double **coord, double center[3] )
{
  int i;
  double  r, temp_r;
  r = 0.0;
  for( i = 0 ; i < N ; i++ ) {
    temp_r = dist_d( coord[i], center, 3 );
    if (temp_r > r) {
      r = temp_r;
    }
  }
  return r;
}

/**
 *
 *
 * @param[in]      *A                3D point.
 * @param[in]      *B                3D point.
 * @param[out]     *C                3D point.
 *
 * @return \f$\angle ABC\f$ in radians
 */
double angle_ABC(double *A,double *B, double *C)
{
	int k;
	double u[3],w[3];
	for(k=0;k<3;k++) { u[k] = A[k] -	B[k]; }
	for(k=0;k<3;k++) { w[k] = C[k] -	B[k]; }
	normalize_d(u,3);
	normalize_d(w,3);
	return acos(scal_d(u,w,3));
}


/**
 *
 *
 * @param[in,out] **coord        N 3D points to be rotated.
 * @param[in]       N            Number of points.
 * @param[in]       center       3D point used as center of rotation.
 * @param[in]       axis         Vector around which the rotation is performed.
 * @param[in]       theta        rotation angle.
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

  if ( coord == NULL || coord[0]==NULL) {
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
  for ( i = N-1; i >=0 ; i-- ) {
    coord[i][0] -=center[0];
    coord[i][1] -=center[1];
    coord[i][2] -=center[2];
  }
  for ( i = N-1; i >=0 ; i-- ) {
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
 *
 * @param[in]      *v            3D vector from which base is calculated.
 *
 * @return Orhonormal basis in 3D with *v as one of coordinates
 */
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
	for(i=0;i<n;i++) {
		vw[0][i]=gsl_vector_get(v,i);
		norm_v2+=vw[0][i]*vw[0][i];
	}
	norm_v=sqrt(norm_v2);
	//construct basis
	for(i=1;i<n;i++) {
		if(i==k) {
			vw[k][0]=1.;
			for(j=1;j<n;j++) { vw[k][j]=0; }
		} else {
			for(j=0;j<n;j++) { vw[i][j]=0.;}
			vw[i][i]=1.;
		}
	}
	//orthonormalize and save basis
	for(k=0;k<n;k++) {
		vw[0][k]/=norm_v;
		gsl_matrix_set(B,0,k,vw[0][k]);
	}
	for(i=0;i<n;i++) {
		if(i>0){
			norm_v2=0;
			for(k=0;k<n;k++) { norm_v2+=vw[i][k]*vw[i][k]; }
			norm_v=sqrt(norm_v2);
			for(k=0;k<n;k++) {
				vw[i][k]/=norm_v;
				gsl_matrix_set(B,i,k,vw[i][k]);
			}
		}
		for(j=i+1;j<n;j++) {
			a=0;
			for(k=0;k<n;k++) { a+=vw[j][k]*vw[i][k]; }
			for(k=0;k<n;k++) { vw[j][k]-=a*vw[i][k]; }
		}
	}
	return B;
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
double dihedralangle_ABCD(const double *atom_1, const double *atom_2, const double *atom_3, const double *atom_4)
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
