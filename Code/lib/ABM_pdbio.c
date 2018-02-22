#define _GNU_SOURCE
#include "ABM_pdbio.h"
#include "my_memory.h"
#include <string.h>
#include <ctype.h>
#include <stdio.h>

const pdb_atom PDBIO_default_atom_entry = { 0, //serial
	"XXXX", 				//name
	'X',						//altLoc
	"XXX",					//resName
	'X',						//chainID
	0,							//resSeq
	'X', 						//iCode
	0.0,0.0,0.0,		//x,y,z
	0.0,0.0,				//occupancy,tempFactor
	"XX",						//element
	"XX"};					//charge


/*
int main (int argc, char **argv)
{
	FILE *in,*out;
	int err;
	double **coord=NULL;
	char **at_names=NULL;
	char **resNames=NULL;
	char *altLocs=NULL;
	int *serials=NULL;
	int *resSeqs=NULL;
	pdb_atom props;
	pdb_atom_list *atomsB=NULL;
	pdb_atom_list *atoms=NULL; 

	in=fopen(argv[1],"r");
	out=fopen(argv[2],"w");
	err=PDBIO_read_all_atoms (in, &atoms);

	pdb_atom_list *new_atom;
	int N;
	err=PDBIO_atomlist_get_coord(atoms,&N,&coord);
	err=PDBIO_atomlist_get_names(atoms,&N,&at_names);
	err=PDBIO_atomlist_get_resNames(atoms,&N,&resNames);
	err=PDBIO_atomlist_get_serials(atoms,&N,&serials);
	err=PDBIO_atomlist_get_resSeq	(atoms,&N,&resSeqs);
	err=PDBIO_atomlist_get_altLoc	(atoms,&N,&altLocs);
	props=atoms->values;
	props.occupancy=1.0;
	props.tempFactor=1.0;
	strncpy(props.element,"  ",2);
	strncpy(props.charge,"  ",2);
	//err=PDBIO_write_all_atoms(out, atoms);
	err=PDBIO_coord2atomlist (coord, serials,at_names, altLocs,resNames, resSeqs, N, props, &atomsB);
	//err=PDBIO_coord2atomlist (coord, NULL,at_names, NULL,resNames, resSeqs, N, props, &atomsB);
	err=PDBIO_write_all_atoms(out, atomsB);

	fclose(in);
  fclose(out);

	return 0;
}
*/


inline void strip_spaces(char *s)
{
	//Quite clever trimming loop I've found at
	//http://stackoverflow.com/questions/7775138/strip-whitespace-from-a-string-in-place
	//The termination condition tests that s[i]==NULL and set s[j]=s[i] in one swoop.
	//Tested OK. LT 01.08.16
	for ( size_t i=0, j=0;(s[j]=s[i]); j+=!isspace(s[i++]));
}


int PDBIO_read_atom ( char * line, pdb_atom *at)
{
	//Tested OK. LT 02.08.16
	char *s;
	char buff_line[81];
	if (line == NULL)
	{
		return PDBIO_ERROR_EMPTY_LINE;
	}
	/*if (strlen(line)<80)
	{
		return PDBIO_ERROR_SHORT_LINE;
	}
	*/
	if ( at == NULL)
	{
		return PDBIO_ERROR_NULL_PTR;
	}
	memset(buff_line,  ' ',80);
	buff_line[80]='\0';
	memset(at->name,   '\0',sizeof(at->name));
	memset(at->resName,'\0',sizeof(at->resName));
	memset(at->element,'\0',sizeof(at->element));
	memset(at->charge, '\0',sizeof(at->charge));
	snprintf(buff_line,80,"%s",line);
	if(buff_line[21]!=' ') //STD behavior: chainID is a letter.
	{
		sscanf(buff_line,"%*6c%5d%*c%4c%c"     //1-17
				"%3c %c%4d%c   "    //18-30 (spaces are important!)
				"%8f%8f%8f"  				//31-54 (coords)
				"%6f%6f%*10c%2c%2c",&at->serial,at->name,&at->altLoc,
				at->resName,&at->chainID,&at->resSeq,&at->iCode,
				&at->x,&at->y,&at->z,
				&at->occupancy,&at->tempFactor,
				at->element,at->charge);
	}
	else
	{
		sscanf(buff_line,"%*6c%5d%*c%4c%c"     //1-17
				"%3c  %4d%c   "    //18-30 (spaces are important!)
				"%8f%8f%8f"  				//31-54 (coords)
				"%6f%6f%*10c%2c%2c",&at->serial,at->name,&at->altLoc,
				at->resName,&at->resSeq,&at->iCode,
				&at->x,&at->y,&at->z,
				&at->occupancy,&at->tempFactor,
				at->element,at->charge);
		at->chainID=' ';
	}
	strip_spaces(at->name);
	strip_spaces(at->resName);
	strip_spaces(at->element);
	strip_spaces(at->charge);
	return PDBIO_SUCCESS;
}

int PDBIO_write_atom( pdb_atom *at, char *line)
{
	//Tested OK. LT 01.08.16
	if (line == NULL)
	{
		return PDBIO_ERROR_EMPTY_LINE;
	}
	//from http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOMge
	//Alignment of one-letter atom name such as C starts at column 14, 
	//while two-letter atom name such as FE starts at column 13.
	if(!isdigit(at->name[0]) && !isalnum(at->name[3]))
	{
		snprintf(line,80,"ATOM  %5d  %-3s%c"    //1-17
				"%3s %c%4d%c   "    		//18-30 (space are important!)
				"%8.3f%8.3f%8.3f"  			//31-54 (coords)
				"%6.2f%6.2f          "  //55-72 (space are important!)
				"%2s%2s",at->serial,at->name,at->altLoc,
				at->resName,at->chainID,at->resSeq,at->iCode,
				at->x,at->y,at->z,
				at->occupancy,at->tempFactor,
				at->element,at->charge);
	}
	else
	{
		snprintf(line,80,"ATOM  %5d %4s%c"      //1-17
				"%3s %c%4d%c   "    		//18-30 (space are important!)
				"%8.3f%8.3f%8.3f"  			//31-54 (coords)
				"%6.2f%6.2f          "  //55-72 (space are important!)
				"%2s%2s",at->serial,at->name,at->altLoc,
				at->resName,at->chainID,at->resSeq,at->iCode,
				at->x,at->y,at->z,
				at->occupancy,at->tempFactor,
				at->element,at->charge);
	}

	return PDBIO_SUCCESS;
}
int PDBIO_read_all_atoms (FILE * pdb_in, pdb_atom_list ** atoms)
{
	//Tested OK. LT 01.08.16 -- NOT VALID ANYMORE
	char *line;
	size_t line_len;
	char line_ID[5];
	int pdberr;
	int cnt;
	pdb_atom atom;
	pdb_atom_list *old_atoms=*atoms;
	pdb_atom_list *tmp_atoms=NULL;
	pdb_atom_list *tmp_atoms_tail=NULL;
	//pdb_atom_list *new_atom;
	cnt=0;
	if(pdb_in==NULL)
	{
		return PDBIO_ERROR_NULL_PTR ;
	}
	line=(char*)malloc(2048*sizeof(char));
	while((getline(&line,&line_len,pdb_in))!=-1)
	{
		strncpy(line_ID,line,4);
		line_ID[4]='\0';
		if (strncmp(line_ID,"ATOM",4)==0)
		{
			pdberr=PDBIO_read_atom(line,&atom);
			PDBIO_atomlist_append(&tmp_atoms_tail,atom);
			if(NULL==tmp_atoms){ tmp_atoms=tmp_atoms_tail; }
			cnt++;

			/*
			if(cnt==0)
			{
				tmp_atoms=(pdb_atom_list*)malloc(sizeof(pdb_atom_list));
				new_atom=tmp_atoms;
			}
			else
			{
				new_atom->next_atom=(pdb_atom_list*)malloc(sizeof(pdb_atom_list));
				new_atom=new_atom->next_atom;
			}
			pdberr=PDBIO_read_atom(line,&new_atom->values);
			if(pdberr !=0 )
			{
				fprintf(stderr,"failed reading file! error code %d\n",pdberr);
				return pdberr;
			}
			new_atom->next_atom=NULL;
			cnt++;
			*/
		}
	}
	if (*atoms==NULL)
	{
		*atoms=tmp_atoms;
	}
	else if ((*atoms)->next_atom != NULL )
	{
		return PDBIO_ERROR_OVERWRITE;
	}
	else
	{
		(*atoms)->next_atom=tmp_atoms;
	}
	free(line);
	return PDBIO_SUCCESS;
}

int PDBIO_write_all_atoms ( FILE *pdb_out, pdb_atom_list *atoms)
{
	//Tested OK. LT 01.08.16
	int pdberr;
	char *line=(char*)malloc(81*sizeof(char));
	pdb_atom_list *new_atom;
	//new_atom=atoms;
	//while( new_atom!=NULL)
	for(new_atom=atoms;new_atom;new_atom=new_atom->next_atom)
	{
		pdberr=PDBIO_write_atom(&new_atom->values,line);
		if (pdberr !=0 )
		{
			return pdberr;
		}
		fprintf(pdb_out,"%s\n",line);
		//new_atom=new_atom->next_atom;
	}
	free(line);
	return PDBIO_SUCCESS;
}

/*
 * push a new element on top of stack
 */
void PDBIO_atomlist_push ( pdb_atom_list **head, pdb_atom el )
{
	//NOT TESTED
	pdb_atom_list * l = ( pdb_atom_list * ) malloc(sizeof(pdb_atom_list));
	if ( l == NULL )
	{
		failed("In Function PDBIO_push. Failed to allocate memory in push.ABORT");
	}
	l->values=el;
	l->prev_atom = NULL;
	l->next_atom = *head;
	if(NULL != *head)
	{
		(*head)->prev_atom=l;
	}
	*head = l;
}

/*
 * recover top element from top of stack
 */
void PDBIO_atomlist_pop ( pdb_atom_list ** head ,pdb_atom *el)
{
	//NOT TESTED
  pdb_atom_list * l;
  if ( *head == NULL )
  {
    fprintf(stderr,"In function PDBIO_atomlist_pop.\n Trying to pop an empty stack.\nABORT\n");
    exit(1);
  }
  l = *head;

  *el = l->values;

  *head = l->next_atom;
	(*head)->prev_atom=NULL;

  free(l);
}

/*
 * append a new element at the end of the list
 */
void PDBIO_atomlist_append ( pdb_atom_list **tail, pdb_atom el )
{
	//NOT TESTED
	pdb_atom_list * l = ( pdb_atom_list * ) malloc(sizeof(pdb_atom_list));
	if ( l == NULL )
	{
		failed("In Function PDBIO_append. Failed to allocate memory in append.ABORT");
	}
	l->values=el;
	l->next_atom=NULL;
	if ( NULL == *tail)
	{
		l->prev_atom=NULL;
		*tail=l;
	}
	else
	{
		(*tail)->next_atom=l;
		l->prev_atom=*tail;
		*tail=l;
	}
}
/*
 * remove an element from the list
 */
void PDBIO_atomlist_remove ( pdb_atom_list ** head ,pdb_atom_list *obj)
{
	pdb_atom_list *pp,*pn;
	if ( *head == obj )
	{
		*head=obj->next_atom;
		(*head)->prev_atom=NULL;
	}
	else
	{
		pp=obj->prev_atom;
		pn=obj->next_atom;
		pp->next_atom=pn;
		if (NULL != pn )
		{
			pn->prev_atom=pp;
		}
		free(obj);
	}
}
int PDBIO_atomlist_get_coord		 ( pdb_atom_list *atoms, int *N_at, double	 	 ***coord)
{
	//Tested OK. LT 01.08.16
	int N;
	int i;
	pdb_atom_list *new_atom;
	double **c;

	if ( *coord != NULL)
	{
		return PDBIO_ERROR_OVERWRITE;
	}
	if (N_at ==NULL || atoms==NULL)
	{
		return PDBIO_ERROR_NULL_PTR;
	}

	N=0;
	for(new_atom=atoms;new_atom;new_atom=new_atom->next_atom) { N++;}
	*N_at=N;
	c=d2t(N,3);
	i=0;
	for(new_atom=atoms;new_atom;new_atom=new_atom->next_atom) 
	{
		c[i][0]=new_atom->values.x;
		c[i][1]=new_atom->values.y;
		c[i][2]=new_atom->values.z;
		i++;
	}
	*coord=c;
	return PDBIO_SUCCESS;
}
int PDBIO_atomlist_get_serials ( pdb_atom_list *atoms, int *N_at, int **serials)
{
	//Tested OK. LT 01.08.16
	int N;
	int i;
	int *s;
	pdb_atom_list *new_atom;

	if ( *serials != NULL)
	{
		return PDBIO_ERROR_OVERWRITE;
	}
	if (N_at ==NULL || atoms==NULL)
	{
		return PDBIO_ERROR_NULL_PTR;
	}
	N=0;
	for(new_atom=atoms;new_atom;new_atom=new_atom->next_atom) { N++;}
	*N_at=N;
	s=i1t(N);
	i=0;
	for(new_atom=atoms;new_atom;new_atom=new_atom->next_atom) 
	{
		s[i]=new_atom->values.serial;
		i++;
	}
	*serials=s;
	return PDBIO_SUCCESS;
}
int PDBIO_atomlist_get_names ( pdb_atom_list *atoms, int *N_at, char ***names)
{
	//Tested OK. LT 01.08.16
	int N;
	int i;
	pdb_atom_list *new_atom;
	char **n;

	if ( *names != NULL)
	{
		return PDBIO_ERROR_OVERWRITE;
	}
	if (N_at ==NULL || atoms==NULL)
	{
		return PDBIO_ERROR_NULL_PTR;
	}
	N=0;
	for(new_atom=atoms;new_atom;new_atom=new_atom->next_atom) { N++;}
	*N_at=N;
	n=c2t(N,5);
	i=0;
	for(new_atom=atoms;new_atom;new_atom=new_atom->next_atom) 
	{
		strncpy(n[i],new_atom->values.name,5);
		i++;
	}
	*names=n;
	return PDBIO_SUCCESS;
}
int PDBIO_atomlist_get_resNames ( pdb_atom_list *atoms, int *N_at, char ***resNames)
{
	//Tested OK. LT 01.08.16
	int N;
	int i;
	char **r;
	pdb_atom_list *new_atom;

	if ( *resNames != NULL)
	{
		return PDBIO_ERROR_OVERWRITE;
	}
	if (N_at ==NULL || atoms==NULL)
	{
		return PDBIO_ERROR_NULL_PTR;
	}
	N=0;
	for(new_atom=atoms;new_atom;new_atom=new_atom->next_atom) { N++;}
	*N_at=N;
	r=c2t(N,4);
	i=0;
	for(new_atom=atoms;new_atom;new_atom=new_atom->next_atom) 
	{
		strncpy(r[i],new_atom->values.resName,4);
		i++;
	}
	*resNames=r;
	return PDBIO_SUCCESS;
}
int PDBIO_atomlist_get_altLoc ( pdb_atom_list *atoms, int *N_at, char **altLocs)
{
	//Tested OK. LT 01.08.16
	int N;
	int i;
	char *a;
	pdb_atom_list *new_atom;

	if ( *altLocs != NULL)
	{
		return PDBIO_ERROR_OVERWRITE;
	}
	if (N_at ==NULL || atoms == NULL)
	{
		return PDBIO_ERROR_NULL_PTR;
	}
	N=0;
	for(new_atom=atoms;new_atom;new_atom=new_atom->next_atom) { N++;}
	*N_at=N;
	a=c1t(N);
	i=0;
	for(new_atom=atoms;new_atom;new_atom=new_atom->next_atom) 
	{
		a[i]=new_atom->values.altLoc;
		i++;
	}
	*altLocs=a;
	return PDBIO_SUCCESS;
}
int PDBIO_atomlist_get_resSeq ( pdb_atom_list *atoms, int *N_at, int **resSeqs)
{
	//Tested OK. LT 01.08.16
	int N;
	int i;
	int *r;
	pdb_atom_list *new_atom;

	if ( *resSeqs != NULL)
	{
		return PDBIO_ERROR_OVERWRITE;
	}
	if (N_at ==NULL || atoms==NULL)
	{
		return PDBIO_ERROR_NULL_PTR;
	}
	N=0;
	for(new_atom=atoms;new_atom;new_atom=new_atom->next_atom) { N++;}
	*N_at=N;
	r=i1t(N);
	i=0;
	for(new_atom=atoms;new_atom;new_atom=new_atom->next_atom) 
	{
		r[i]=new_atom->values.resSeq;
		i++;
	}
	*resSeqs=r;
	return PDBIO_SUCCESS;
}
int PDBIO_atomlist_get_occupancy ( pdb_atom_list *atoms, int *N_at, double		**occupancies)
{
	int N;
	int i;
	pdb_atom_list *new_atom;
	double *o;

	if ( *occupancies != NULL)
	{
		return PDBIO_ERROR_OVERWRITE;
	}
	if (N_at ==NULL || atoms==NULL)
	{
		return PDBIO_ERROR_NULL_PTR;
	}
	N=0;
	for(new_atom=atoms;new_atom;new_atom=new_atom->next_atom) { N++;}
	*N_at=N;
	o=d1t(N);
	i=0;
	for(new_atom=atoms;new_atom;new_atom=new_atom->next_atom) 
	{
		o[i]=new_atom->values.occupancy;
		i++;
	}
	*occupancies=o;
	return PDBIO_SUCCESS;
}
int PDBIO_atomlist_get_tempFactor ( pdb_atom_list *atoms, int *N_at, double		**tempFactors)
{
	int N;
	int i;
	double *tF;
	pdb_atom_list *new_atom;

	if ( *tempFactors!= NULL)
	{
		return PDBIO_ERROR_OVERWRITE;
	}
	if (N_at ==NULL || atoms==NULL)
	{
		return PDBIO_ERROR_NULL_PTR;
	}
	N=0;
	for(new_atom=atoms;new_atom;new_atom=new_atom->next_atom) { N++;}
	*N_at=N;
	tF=d1t(N);
	i=0;
	for(new_atom=atoms;new_atom;new_atom=new_atom->next_atom) 
	{
		tF[i]=new_atom->values.occupancy;
		i++;
	}
	*tempFactors=tF;
	return PDBIO_SUCCESS;
}

/*int PDBIO_atomlist2coord (  pdb_atom_list *atoms, pdb_atom *chain_prop,int **serial,char ***name, char **altLoc, char ***resName, int **resSeq, int *N_at, double ***coord)
{
	pdb_atom_list *new_atom;
	int N;
	for(new_atom=atoms;new_atom;new_atom=new_atom->next_atom) { N++;}
	coord=d2t(N,3);
	at_names=c2t(N,5);
	resNames=c2t(N,4);
	serials=i1t(N);
	resSeqs=i1t(N);
	altLocs=c1t(N);
	int i=0;
	for(new_atom=atoms;new_atom;new_atom=new_atom->next_atom) 
	{
		coord[i][0]=new_atom->values.x;
		coord[i][1]=new_atom->values.y;
		coord[i][2]=new_atom->values.z;
		strncpy(at_names[i],new_atom->values.name,5);
		strncpy(resNames[i],new_atom->values.resName,4);
		serials[i]=new_atom->values.serial;
		altLocs[i]=new_atom->values.altLoc;
		resSeqs[i]=new_atom->values.resSeq;
		i++;
	}
}*/
int PDBIO_coord2atomlist (double **coord, int *serial,char **at_name, char *altLoc, char **resName, int *resSeq, int N_at, pdb_atom chain_prop, pdb_atom_list **atoms)
{
	//Tested OK. LT 01.08.16
	int i;
	pdb_atom atom=chain_prop;
	if ( coord == NULL)
	{
		return PDBIO_ERROR_NULL_PTR;
	}
	//list atoms is treated as a stack: first in is last out.
	for(i=N_at-1; i>=0; i--)
	{
		atom.x=coord[i][0];
		atom.y=coord[i][1];
		atom.z=coord[i][2];
		if (altLoc!=NULL) {atom.altLoc=altLoc[i];}
		else {atom.altLoc=' ';}
		if (serial!=NULL){ atom.serial=serial[i];}
		else {atom.serial=i+1;}
		atom.resSeq=resSeq[i];
		strncpy(atom.name,at_name[i],5);
		strncpy(atom.resName,resName[i],4);
		snprintf(atom.element,3,"  ");
		atom.element[1]=atom.name[0];
		PDBIO_atomlist_push(atoms,atom);
	}
	return PDBIO_SUCCESS;
}

