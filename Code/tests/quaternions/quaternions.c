#include <stdio.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>

#include "../../lib/quaternions.h"
#include "../../lib/Caterpillar.h"
#include "../../lib/Caterpillar_IO.h"
#include "../../lib/my_geom.h"

void build_peptide_x ( gsl_matrix *pep);
int add_peptide ( cat_prot *p, int I, double phi, double alpha, double psi );
int prot_rescale ( cat_prot *p);
void simple_PDB_write_atom (int serial, char *atName, char *resName, int resSeq, double x, double y, double z, char *line);
void print_peptide_PDB ( FILE * pdbout, gsl_matrix * pep);
int rototran_peptide ( gsl_matrix *pep_r, double phi, double alpha, double psi, gsl_vector *axPhi, gsl_vector *t, gsl_matrix *pep);
int main (int argc, char **argv)
{
	int N_prot;
	cat_prot *p;
	char  ATnames_to_read[5][5]={"N","CA","C","O","H"};
	CATIO_pdb2cat (&p,&N_prot,5,ATnames_to_read,'A',argv[1]);
	CAT_insert_hydrogens(p);
	CATIO_cat2pdb("prova_rescale.pdb","w","",p,1);
	prot_rescale(p);
	CATIO_cat2pdb("prova_rescale.pdb","a","",p,1);
	CAT_prot_free(p);
	return 0;
}

/*
void build_peptide_x ( gsl_matrix *pep)
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
	theta=CAT_angle_NCaC-M_PI_2;
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
	theta= CAT_angle_NCaC-M_PI_2+beta; //mind the sign
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
*/
void build_peptide_x ( gsl_matrix *pep)
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

/*
int rototran_peptide ( gsl_matrix *pep_r, double phi, double alpha, double psi, gsl_vector *axPhi, gsl_vector *t, gsl_matrix *pep)
{
	double q1[4]={0,0,0,0};
	double q2[4]={0,0,0,0};
	double q3[4]={0,0,0,0};
	double qmul[4]={0,0,0,0};
	gsl_vector_view q1_v 	=gsl_vector_view_array(q1,4);
	gsl_vector_view q2_v 	=gsl_vector_view_array(q2,4);
	gsl_vector_view q3_v	=gsl_vector_view_array(q3,4);
	gsl_vector_view qmul_v	=gsl_vector_view_array(qmul,4);
	double axAlpha[3]={0,0,1};
	double axPsi[3];
	gsl_vector_view axAlpha_v =gsl_vector_view_array(axAlpha,3);
	gsl_vector_view axPsi_v		=gsl_vector_view_array(axPsi,3);
	gsl_vector_view orig =gsl_matrix_column(pep,1);
	//Phi
	quaternion_build(&q1_v.vector,axPhi,phi);
	//Alpha
	quaternion_build(&q2_v.vector,&axAlpha_v.vector,alpha);
	//Psi
	gsl_vector_memcpy(&axPsi_v.vector,&orig.vector);
	orig =gsl_matrix_column(pep,0);
	gsl_vector_sub(&axPsi_v.vector, &orig.vector);
	quaternion_build(&q3_v.vector,&axPsi_v.vector,psi);
	//multiply quaternions
	quaternion_mult (&qmul_v.vector,&q2_v.vector,&q1_v.vector);
	gsl_vector_memcpy(&q2_v.vector,&qmul_v.vector);
	quaternion_mult (&qmul_v.vector,&q3_v.vector,&q2_v.vector);
	//build rototranslation
	double rt[16];
	gsl_matrix_view rt_v = gsl_matrix_view_array(rt,4,4);
	rototransl3D_build(&rt_v.matrix,&qmul_v.vector,t);
	//adapt peptide coord for rototransl.
	double gen_coord[4*6];
	double gen_coord_rot[4*6];
	gsl_matrix_view gc_v = gsl_matrix_view_array(gen_coord,4,6);
	gsl_matrix_view gcr_v = gsl_matrix_view_array(gen_coord_rot,4,6);
	gsl_matrix_set_zero(&gc_v.matrix);
	gsl_matrix_view c_v= gsl_matrix_submatrix(&gc_v.matrix,0,0,3,6);
	gsl_matrix_memcpy(&c_v.matrix,pep);
	for(int j=0;j<6;j++){
		gsl_matrix_set(&gc_v.matrix,3,j,1);
	}
	//rototranslate pep
	printf("-----\n");
	for(int i=0;i<4;i++) {
		for(int j=0;j<4;j++) {
			fprintf(stdout, "%6.4f ",gsl_matrix_get(&rt_v.matrix,i,j));
		}
		fprintf(stdout,"\n");
	}
	printf("-----\n");
	for(int i=0;i<4;i++) {
		for(int j=0;j<6;j++) {
			fprintf(stdout, "%6.4f ",gsl_matrix_get(&gc_v.matrix,i,j));
		}
		fprintf(stdout,"\n");
	}
	printf("-----\n");
	int error=gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.,&rt_v.matrix,&gc_v.matrix,0.,&gcr_v.matrix);
	for(int i=0;i<4;i++) {
		for(int j=0;j<6;j++) {
			fprintf(stdout, "%6.4f ",gsl_matrix_get(&gcr_v.matrix,i,j));
		}
		fprintf(stdout,"\n");
	}
	//get coordinates
	c_v= gsl_matrix_submatrix(&gcr_v.matrix,0,0,3,6);
	gsl_matrix_memcpy(pep_r,&c_v.matrix);
	return error;
}
*/

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
		printf("%d alpha: %lf\n",i,alpha[i]/M_PI*180);
		error=add_peptide(p,i,  p->phi[i],M_PI-alpha[i],p->psi[i]);
	}
	compute_dihedrals(p);

	return error;
}


int add_peptide ( cat_prot *p, int I, double phi, double alpha, double psi )
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
	build_peptide_x (&pep_v.matrix); //returns CA-C-O-N-H-CA(I+1)
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


void simple_PDB_write_atom (int serial, char *atName, char *resName, int resSeq, double x, double y, double z, char *line)
{
	snprintf(line,80,"ATOM  %5d  %-3s%c"    //1-17
			"%3s %c%4d%c   "    		//18-30 (space are important!)
			"%8.3f%8.3f%8.3f"  			//31-54 (coords)
			"%6.2f%6.2f          "  //55-72 (space are important!)
			"%2s%2s",serial,atName,' ',
			resName,'A',resSeq,' ',
			x,y,z,
			1.,1.,
			atName," ");

}

void print_peptide_PDB ( FILE * pdbout, gsl_matrix * pep)
{
	char line[90];
	simple_PDB_write_atom (1, "CA", "ALA",1,
			gsl_matrix_get(pep,0,0), gsl_matrix_get(pep,1,0),gsl_matrix_get(pep,2,0),line);
	fprintf(pdbout,"%s\n",line);
	simple_PDB_write_atom (2, "C", "ALA",1,
			gsl_matrix_get(pep,0,1), gsl_matrix_get(pep,1,1),gsl_matrix_get(pep,2,1),line);
	fprintf(pdbout,"%s\n",line);
	simple_PDB_write_atom (3, "O", "ALA",1,
			gsl_matrix_get(pep,0,2), gsl_matrix_get(pep,1,2),gsl_matrix_get(pep,2,2),line);
	fprintf(pdbout,"%s\n",line);
	simple_PDB_write_atom (4, "N", "ALA",1,
			gsl_matrix_get(pep,0,3), gsl_matrix_get(pep,1,3),gsl_matrix_get(pep,2,3),line);
	fprintf(pdbout,"%s\n",line);
	simple_PDB_write_atom (5, "H", "ALA",1,
			gsl_matrix_get(pep,0,4), gsl_matrix_get(pep,1,4),gsl_matrix_get(pep,2,4),line);
	fprintf(pdbout,"%s\n",line);
	simple_PDB_write_atom (6, "CA", "ALA",1,
			gsl_matrix_get(pep,0,5), gsl_matrix_get(pep,1,5),gsl_matrix_get(pep,2,5),line);
	fprintf(pdbout,"%s\n",line);
	fprintf(pdbout,"ENDMDL\n");
	return;
}
