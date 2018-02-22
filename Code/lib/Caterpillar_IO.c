#include "Caterpillar_IO.h"
#include "messages.h"
#include "my_memory.h"

void CATIO_fastadecoder (int *Dec, char *Enc,int Seq_Length)
{
	int i;
	int err5=1;
	for(i=0;i<Seq_Length;i++){
		switch (Enc[i]) {
			case 'A':
				Dec[i]=1;
				break;
			case 'C':
				Dec[i]=2;
				break;
			case 'D':
				Dec[i]=3;
				break;
			case 'E':
				Dec[i]=4;
				break;
			case 'F':
				Dec[i]=5;
				break;
			case 'G':
				Dec[i]=6;
				break;
			case 'H':
				Dec[i]=7;
				break;
			case 'I':
				Dec[i]=8;
				break;
			case 'K':
				Dec[i]=9;
				break;
			case 'L':
				Dec[i]=10;
				break;
			case 'M':
				Dec[i]=11;
				break;
			case 'N':
				Dec[i]=12;
				break;
			case 'P':
				Dec[i]=13;
				break;
			case 'Q':
				Dec[i]=14;
				break;
			case 'R':
				Dec[i]=15;
				break;
			case 'S':
				Dec[i]=16;
				break;
			case 'T':
				Dec[i]=17;
				break;
			case 'V':
				Dec[i]=18;
				break;
			case 'W':
				Dec[i]=19;
				break;
			case 'Y':
				Dec[i]=20;
				break;
			case 'X':
				Dec[i]=1;
				break;
			default:
				fprintf(stderr,"%s:%d Errore nel codice della Senquenza Residuo Enc[%d]=%c unknown\n",__FILE__,__LINE__,i,Enc[i]);
				fflush(stderr);
				failed("error in CATIO_fastadecoder.");
		}
	}
}

void CATIO_cat2pdb 		(  char *filename, char *writemode, char *remark, cat_prot *prot,  int N_prot)
{
	//Tested OK for 1 protein and writemode "w" only. LT 02.08.16
	int i,j,k,l;
	int N_at;
	int N_res;
	int N_atom_per_res;
	char **at_names;
	char **resNames;
	int *resSeqs;
	char CAT_AT[6][6]=CAT_ATnames;
	char CAT_AA[CAT_S][4]=CAT_AADefs;
	int  CAT_AC[CAT_S]=CAT_AACODE;
	char C_IDs[25]="ABCDEFGHIJKLMNOPQRSTUVXYZ";
	char error[1024];
	char start_string[128];
	char end_string[128];
	FILE *pdb_out;
	pdb_atom chain_prop;
	pdb_atom_list *all_atoms=NULL;

	for(i=0;i<N_prot;i++)
	{
		N_at=prot[i].n_atoms;
		N_res=prot[i].n_res;
		N_atom_per_res=prot[i].n_atom_per_res;
		resSeqs=i1t(N_at);
		at_names=c2t(N_at,5);
		resNames=c2t(N_at,4);
		for(j=0;j<N_res;j++)
		{
			l=prot[i].residues[j];
			for(k=0;k<N_atom_per_res;k++)
			{
				resSeqs[j*N_atom_per_res+k]=j;
				snprintf(at_names[j*N_atom_per_res+k],3,"%s",CAT_AT[k]);
				snprintf(resNames[j*N_atom_per_res+k],4,"%s",CAT_AA[CAT_AC[l]]);//asdas
				strip_spaces(at_names[j*N_atom_per_res+k]);
				strip_spaces(resNames[j*N_atom_per_res+k]);
			}
		}
		chain_prop.chainID=C_IDs[i];
		chain_prop.occupancy	=1.0;
		chain_prop.tempFactor	=1.0;
		chain_prop.iCode=' ';
		memset(chain_prop.charge,'\0',sizeof(chain_prop.charge));
		strncpy(chain_prop.charge,"  ",2);
		PDBIO_coord2atomlist( prot[i].coord,NULL, at_names,NULL, resNames,resSeqs, prot[i].n_atoms, chain_prop, &all_atoms);
		free_i1t(resSeqs);
		free_c2t(at_names);
		free_c2t(resNames);
	}
	//OUTPUT on file
	if(writemode[0]=='w')
	{
		start_string[0]='\0';
		sprintf(end_string,"END");
	}
	else if(writemode[0]=='a')
	{
		sprintf(start_string,"MODEL");
		sprintf(end_string,"ENDMDL");
	}
	else
	{
		sprintf(error,"writemode not recognized: %s\n",writemode);
		failed(error);
	}
  if ((pdb_out=fopen(filename,writemode))==NULL)
  {
    sprintf(error,"Could not open file %s (%s).\n",filename,writemode);
    failed(error);
  }
	//
	if(strlen(remark)>0)
	{
		fprintf(pdb_out,"REMARK %s\n",remark);
	}
	//fprintf(pdb_out,"%s\n",start_string);
	PDBIO_write_all_atoms(pdb_out,all_atoms); //should specify writemode!!!
	fprintf(pdb_out,"%s\n",end_string);
	//close
	fclose(pdb_out);
	//free
	pdb_atom_list *pn,*pp=all_atoms;
	while(pp)
	{
		pn=pp->next_atom;
		free(pp);
		pp=pn;
	}
}

void CATIO_pdb2cat  	( cat_prot **prot, int *N_prot, int N_atom_types, char ATnames[N_atom_types][5], char Loc_keep, char *filename)
{
	//Tested OK 12/07/17
	//Note: does not insert hydrogens, does not rescale, does not compute dihedrals
	pdb_atom_list *atoms=NULL;
	pdb_atom_list *pn,*pp;
	char 	error[1024];
	int 	IO_err;
	FILE *pdb_in;
	int  	N_at;
	int delete;
	char *s;
	double **buff_coord=NULL;
	//char   **buff_names;
	//char   **buff_resNames;
	char  AA_3lett_code[CAT_S][4]=CAT_AADefs;
	int 	CAT_AAcode[CAT_S]=CAT_AACODE;
	char  cat_ATnames[NATOM][6]=CAT_ATnames;
	int flag_H_insert=1; //default: insert hydrogens
	cat_prot *p;
	for(int i=0;i<NATOM;i++) {strip_spaces(ATnames[i]);}
	if ( *prot !=NULL)
	{
		sprintf(error,"Possible overwrite of an already allocated cat_prot pointer. In function CATIO_pdb2cat\n");
	}
  if ((pdb_in=fopen(filename,"r"))==NULL)
	{
		sprintf(error,"In function CATIO_pdb2cat. Can not open file %s for reading.",filename);
		failed(error);
	}
	IO_err=PDBIO_read_all_atoms ( pdb_in, &atoms);
	pn=atoms;
	while(pn)
	{
		delete=1;
		strip_spaces(pn->values.name);
		//keep only the atom types given in input to the function
		for(int i=0;i<N_atom_types;i++)
		{
			if(strncmp(pn->values.name,ATnames[i],5)==0)
			{
				delete=0;
				break;
			}
		}
		//keep only the altLoc given in input (if specified in the file)
		if(pn->values.altLoc!=' ' && pn->values.altLoc!=Loc_keep)
		{
			delete=1;
		}
		if(delete) 
		{
			PDBIO_atomlist_remove ( &atoms , pn);
		}
		pn=pn->next_atom;
	}
	//count how many residues we have (pdb may not have hydrogens).
	int n_res=0;
	for(pn=atoms;pn;pn=pn->next_atom)
	{
		if(strncmp(pn->values.name,cat_ATnames[ATOM_CA],5)==0)
		{
			n_res++;
		}
	}
	IO_err=PDBIO_atomlist_get_coord (atoms,	&N_at,&buff_coord);
	if(N_atom_types==6)
	{
 		p=CAT_prot_alloc (n_res, 6);
	}
	else
	{
 		p=CAT_prot_alloc (n_res, 5);
	}
	int cnt[N_atom_types];
	for(int i=0;i<N_atom_types;i++){cnt[i]=0;}
	for(pn=atoms;pn;pn=pn->next_atom)
	{
		for(int i=0;i<N_atom_types;i++)
		{
			if(strncmp(pn->values.name,ATnames[i],5)==0)
			{
				if(strncmp(pn->values.name,cat_ATnames[ATOM_N],5)==0)
				{
					p->N[cnt[i]][0]=pn->values.x;
					p->N[cnt[i]][1]=pn->values.y;
					p->N[cnt[i]][2]=pn->values.z;
				}
				if(strncmp(pn->values.name,cat_ATnames[ATOM_C],5)==0)
				{
					p->C[cnt[i]][0]=pn->values.x;
					p->C[cnt[i]][1]=pn->values.y;
					p->C[cnt[i]][2]=pn->values.z;
				}
				if(strncmp(pn->values.name,cat_ATnames[ATOM_O],5)==0)
				{
					p->O[cnt[i]][0]=pn->values.x;
					p->O[cnt[i]][1]=pn->values.y;
					p->O[cnt[i]][2]=pn->values.z;
				}
				if(strncmp(pn->values.name,ATnames[ATOM_H],5)==0)
				{
					p->H[cnt[i]][0]=pn->values.x;
					p->H[cnt[i]][1]=pn->values.y;
					p->H[cnt[i]][2]=pn->values.z;
					flag_H_insert=0;
				}
				if(strncmp(pn->values.name,cat_ATnames[ATOM_CA],5)==0)
				{
					p->CA[cnt[i]][0]=pn->values.x;
					p->CA[cnt[i]][1]=pn->values.y;
					p->CA[cnt[i]][2]=pn->values.z;
					//set the residue too.
					int fail=1;
					for(int j=0;j<CAT_S;j++)
					{
						if(strncmp(pn->values.resName,AA_3lett_code[j],3)==0)
						{
							fail=0;
							p->residues[cnt[i]]=CAT_AAcode[j];
							break;
						}
					}
					if( fail )
					{
						sprintf(error,"Error: residue %s can not be identified!\n",pn->values.resName);
						failed(error);
					}
				}
				if(p->CB!=NULL && strncmp(pn->values.name,cat_ATnames[ATOM_CB],5)==0)
				{
					p->CB[cnt[i]][0]=pn->values.x;
					p->CB[cnt[i]][1]=pn->values.y;
					p->CB[cnt[i]][2]=pn->values.z;
				}
				cnt[i]++;
			}
		}
	}
	for(int i=0;i<N_atom_types;i++)
	{
		if(cnt[i]!=n_res)
		{
			sprintf(error,"Error while counting atoms!! Atom %d appears in more than one residue. Atoms:%d n_res: %d",i,cnt[i],n_res);
			failed(error);
		}
	}

	*prot=p;
	free_d2t(buff_coord);
}

void CATIO_pdb2cat_keep_CB  	( cat_prot **prot, int *N_prot, char *filename)
{
	pdb_atom_list *atoms=NULL;
	pdb_atom_list *pn,*pp;
	char 	error[1024];
	int 	IO_err;
	FILE *pdb_in;
	size_t n_atom_per_res=NATOM;
	int  	N_at;
	int delete;
	char *s;
	double **buff_coord=NULL;
	//char   **buff_names;
	//char   **buff_resNames;
	char  AA_3lett_code[CAT_S][4]=CAT_AADefs;
	char  FASTA[CAT_S]=CAT_FASTA;
	int 	CAT_AAcode[CAT_S]=CAT_AACODE;
	char  ATnames[NATOM][5]=CAT_ATnames;
	cat_prot *p;
	for(int i=0;i<NATOM;i++) {strip_spaces(ATnames[i]);}
	if ( *prot !=NULL)
	{
		sprintf(error,"Possible overwrite of an already allocated cat_prot pointer. In function CATIO_pdb2cat\n");
	}
  if ((pdb_in=fopen(filename,"r"))==NULL)
	{
		sprintf(error,"In function CATIO_pdb2cat. Can not open file %s for reading.",filename);
		failed(error);
	}
	IO_err=PDBIO_read_all_atoms ( pdb_in, &atoms);
	pn=atoms;
	pp=NULL;
	while(pn)
	{
		delete=1;
		strip_spaces(pn->values.name);
		for(int i=0;i<NATOM;i++)
		{
			if(strncmp(pn->values.name,ATnames[i],5)==0)
			{
				delete=0;
			}
		}
		if(delete)
		{
			//I should write a simple function for this.
			if(pn==atoms)
			{
				atoms=pn->next_atom;
				free(pn);
				pn=atoms;
			}
			else
			{
				pp->next_atom=pn->next_atom;
				free(pn);
				pn=pp->next_atom;
			}
		}
		else
	 	{
			pp=pn;
			pn=pn->next_atom;
		}
	}
	//count how many residues we have (pdb may not have hydrogens).
	int n_res=0;
	for(pn=atoms;pn;pn=pn->next_atom)
	{
		if(strncmp(pn->values.name,ATnames[ATOM_CA],5)==0)
		{
			n_res++;
		}
	}
	IO_err=PDBIO_atomlist_get_coord (atoms,	&N_at,&buff_coord);
 	p=CAT_prot_alloc (n_res, n_atom_per_res);
	int cnt[NATOM];
	for(int i=0;i<NATOM;i++){cnt[i]=0;}
	for(pn=atoms;pn;pn=pn->next_atom)
	{
		//s=pn->values.name;
		if(strncmp(pn->values.name,ATnames[ATOM_N],5)==0)
		{
			p->N[cnt[ATOM_N]][0]=pn->values.x;
			p->N[cnt[ATOM_N]][1]=pn->values.y;
			p->N[cnt[ATOM_N]][2]=pn->values.z;
			cnt[ATOM_N]++;
		}
		if(strncmp(pn->values.name,ATnames[ATOM_C],5)==0)
		{
			p->C[cnt[ATOM_C]][0]=pn->values.x;
			p->C[cnt[ATOM_C]][1]=pn->values.y;
			p->C[cnt[ATOM_C]][2]=pn->values.z;
			cnt[ATOM_C]++;
		}
		if(strncmp(pn->values.name,ATnames[ATOM_O],5)==0)
		{
			p->O[cnt[ATOM_O]][0]=pn->values.x;
			p->O[cnt[ATOM_O]][1]=pn->values.y;
			p->O[cnt[ATOM_O]][2]=pn->values.z;
			cnt[ATOM_O]++;
		}
		/*
		if(strncmp(pn->values.name,ATnames[ATOM_H],5)==0)
		{
			p->H[cnt[ATOM_H]][0]=pn->values.x;
			p->H[cnt[ATOM_H]][1]=pn->values.y;
			p->H[cnt[ATOM_H]][2]=pn->values.z;
			cnt[ATOM_H]++;
		}
		*/
		//Insert the Hydrogens!
		if(strncmp(pn->values.name,ATnames[ATOM_CA],5)==0)
		{
			p->CA[cnt[ATOM_CA]][0]=pn->values.x;
			p->CA[cnt[ATOM_CA]][1]=pn->values.y;
			p->CA[cnt[ATOM_CA]][2]=pn->values.z;
			//set the residue too.
			int fail=1;
			for(int i=0;i<CAT_S;i++)
			{
				if(strncmp(pn->values.resName,AA_3lett_code[i],3)==0)
				{
					fail=0;
					p->residues[cnt[ATOM_CA]]=CAT_AAcode[i];
					break;
				}
			}
			if( fail )
			{
				sprintf(error,"Error: residue %s can not be identified!\n",pn->values.resName);
				failed(error);
			}
			cnt[ATOM_CA]++;
		}
		if(p->CB!=NULL && strncmp(pn->values.name,ATnames[ATOM_CB],5)==0)
		{
			p->CB[cnt[ATOM_CB]][0]=pn->values.x;
			p->CB[cnt[ATOM_CB]][1]=pn->values.y;
			p->CB[cnt[ATOM_CB]][2]=pn->values.z;
			cnt[ATOM_CB]++;
		}
	}
	for(int i=0;i<NATOM;i++)
	{
		if(cnt[i]!=n_res && i!=ATOM_H)
		{
			sprintf(error,"Error while counting atoms!! Atom %d appears in more than one residue. Atoms:%d n_res: %d",i,cnt[i],n_res);
			failed(error);
		}
	}
	//inser Hydrogens and rescale
	CAT_rescale(p);
	CAT_insert_hydrogens (p);
	CAT_rescale(p);

	*prot=p;
	free_d2t(buff_coord);
}


void CATIO_cat2oldbin 			( char *filename, char *writemode, cat_prot *prot, int N_prot )
{
	double x,y,z;
	int i,j,pr,id;
	FILE *bin_out=NULL;
	char error[1024];
	double **atoms[NATOM];

  if ((bin_out=fopen(filename,writemode))==NULL)
	{
		sprintf(error,"In function CATIO_cat2bin. Can not open file %s for writing with write mode %s",filename, writemode);
		failed(error);
	}
	for(pr=0;pr<N_prot;pr++)
	{
		atoms[ATOM_N]	=prot[pr].N;
		atoms[ATOM_CA]=prot[pr].CA;
		atoms[ATOM_C]	=prot[pr].C;
		atoms[ATOM_O]	=prot[pr].O;
		atoms[ATOM_H]	=prot[pr].H;
		if(prot[pr].CB!=NULL)
		{
			atoms[ATOM_CB]=prot[pr].CB;
		}
		for(i=0;i<prot[pr].n_res;i++)
		{
			for(j=ATOM_N;j<NATOM;j++)
			{
				x=atoms[j][i][0];
				y=atoms[j][i][1];
				z=atoms[j][i][2];
				id=j;
				fwrite(&x,sizeof(double),	1,bin_out);
				fwrite(&y,sizeof(double),	1,bin_out);
				fwrite(&z,sizeof(double),	1,bin_out);
				fwrite(&id,sizeof(int),		1,bin_out);
			}
		}
	}
	fclose(bin_out);
}
void CATIO_oldbin2cat  		( cat_prot *prot, int N_prot, char *filename)
{
	//Compatibility version. N_prot is passed to the function. It would make
	//much more sense to save it on the binary file. Also prot must already be
	//allocate with a valid number of atoms and residues, which will be read from
	//the file.
	double x=0.,y=0.,z=0.0;
	int i=0,pr=0,id=0,sorted_index=0;
	int res=0;
	FILE *bin_in=NULL;
	char error[1024];
 if ((bin_in=fopen(filename,"r"))==NULL)
	{
		sprintf(error,"In function CATIO_bin2cat. Can not open file %s for reading",filename);
		failed(error);
	}

	for(pr=0;pr<N_prot;pr++)
	{
		for(i=0;i<prot[pr].n_atoms;i++)
		{
			fread(&x,sizeof(double),1,bin_in);
			fread(&y,sizeof(double),1,bin_in);
			fread(&z,sizeof(double),1,bin_in);
			fread(&id,sizeof(int),	1,bin_in);
			res=(int)(i/prot[pr].n_atom_per_res)*prot[pr].n_atom_per_res;
			switch(id){
				case ATOM_N:
					sorted_index=res+ATOM_N; // Place the ATOMs in the right order
					break;
				case ATOM_CA:
					sorted_index=res+ATOM_CA;
					break;
				case ATOM_C:
					sorted_index=res+ATOM_C;
					break;
				case ATOM_O:
					sorted_index=res+ATOM_O;
					break;
				case ATOM_H:
					sorted_index=res+ATOM_H;
					break;
				case ATOM_CB:
					sorted_index=res+ATOM_CB;
					break;
			}
		 prot[pr].coord[sorted_index][0]=x;
		 prot[pr].coord[sorted_index][1]=y;
		 prot[pr].coord[sorted_index][2]=z;
		}
	}
	fclose(bin_in);
}

