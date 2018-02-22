#ifndef __MY_QUATERNIONS_HDR
#define __MY_QUATERNIONS_HDR
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
int quaternion_build 		(gsl_vector * q, const gsl_vector *v, const double theta );
int quaternion_mult 		( gsl_vector *q21, const gsl_vector *q2, const gsl_vector *q1);
int rototransl3D_build 	( gsl_matrix * rt, const gsl_vector *quat, const gsl_vector * transl);
int rotation3D_build 		( gsl_matrix * rot, const gsl_vector *quat);
int rototrans3D_mult 		( gsl_matrix *rt21, const gsl_matrix *rt2, const gsl_matrix *rt1);
int rotate3D_vector_m			( gsl_vector *rV, const gsl_matrix *r, const gsl_vector *v);
int rototranslate3D_vector_m (gsl_vector *rtV, const gsl_matrix *rt, const gsl_vector *V);
#endif //__MY_QUATERNIONS_HDR
