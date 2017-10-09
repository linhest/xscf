#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>

#include "molecule.h"
#include "basis_set.h"
extern "C" {
#include <cint.h>
}
#define ENV_SIZE 10000
#define MAXLINELENGTH 2000


basis_set::basis_set(molecule & mol){
  char * filename=mol.get_basis_set();
  cout << "Reading basis:" << mol.get_basis_set() << endl;

  int off = PTR_ENV_START; // = 20
  int n=0;
  for ( int i=0;i<mol.natm;i++){
    read_basis_set_atom(filename, mol.atm_label[i], NULL, NULL, i, &off, &n);
    cout << mol.atm_label[i] << " "<< n << endl;
  }
  nbas=n;
  natm=mol.natm;
  atm = new int[mol.natm* ATM_SLOTS];
  bas = new int[nbas * BAS_SLOTS];
  env = new double[ENV_SIZE];
  off = PTR_ENV_START; // = 20
  n=0;
  for ( int i=0;i<mol.natm;i++){
    atm[CHARGE_OF + ATM_SLOTS * i] = mol.Z[i];
    atm[PTR_COORD + ATM_SLOTS * i] = off;
    env[off + 0] =  mol.geom[i*3+0]; 
    env[off + 1] =  mol.geom[i*3+1]; 
    env[off + 2] =  mol.geom[i*3+2]; 
    off+=3;
  }
  for ( int i=0;i<mol.natm;i++){
    read_basis_set_atom(filename, mol.atm_label[i], env, bas, i, &off, &n);
  }  
}

// reads in basis set data
// from file filename
// for atom atm_label with index atom_of
// env,bas,off,n are the libcint arrays
int basis_set::read_basis_set_atom(char * filename, char * atm_label, double * env, int * bas, int atom_of, int *off, int *n){
  ifstream fin;
  bool found=false;
  int nshell_atm=0;
  fin.open(filename);
  if(!fin.good()) 
    throw std::invalid_argument( strcat(filename," not found ") );
  while(!fin.eof()){
    char alabel[4];
    char slabel[2];
    int nexp, nshell;
    double exp, norm, coef, pcoef;
    int ai;
    int num;
    char buf[MAXLINELENGTH];
    fin.getline(buf,MAXLINELENGTH);
    if(strncmp("****",buf,4)==0 ){
      // new entry starts here
      fin.getline(buf,MAXLINELENGTH);
      if(fin.eof())break;
      num=sscanf(buf,"%s %d",alabel,&ai);
      if(num!=2) throw std::invalid_argument( strcat(buf,"alabel int") );
    }else{
      num=sscanf(buf,"%s %d %lf",slabel,&nexp,&norm);
      if(num!=3) throw std::invalid_argument( strcat(buf,"sym nexp norm") );
      if(strcasecmp(alabel,atm_label)==0 ){
	found=true;
	if(strlen(slabel)==1){
	  nshell_atm++;
	  if(bas!=NULL && env!=NULL){
	    bas[ATOM_OF  + BAS_SLOTS * *n]  = atom_of;
	    bas[ANG_OF   + BAS_SLOTS * *n]  = symlabel2ang(slabel[0]);
	    bas[NPRIM_OF + BAS_SLOTS * *n]  = nexp;
	    bas[NCTR_OF  + BAS_SLOTS * *n]  = 1;
	    bas[PTR_EXP  + BAS_SLOTS * *n]  = *off;
	    bas[PTR_COEFF+ BAS_SLOTS * *n]  = *off+nexp;
	  }
	  for(int i=0;i<nexp;i++){
	    fin.getline(buf,MAXLINELENGTH);
	    num=sscanf(buf,"%lf %lf ",&exp,&coef);
	    if(num!=2) throw std::invalid_argument( strcat(buf,"exp coef") );
	    if(bas!=NULL && env!=NULL){
	      env[*off ] = exp;
	      env[*off+nexp]=coef*CINTgto_norm(bas[ANG_OF+BAS_SLOTS*(*n)],exp);
	    }
	    (*off)++;
	  }
	  (*off)+=nexp;
	  (*n)++;
	}else if(strlen(slabel)==2 && strcmp(slabel,"SP")==0){
	  nshell_atm+=2;		
	  if(bas!=NULL && env!=NULL){
	    bas[ATOM_OF  + BAS_SLOTS * *n]  = atom_of;
	    bas[ANG_OF   + BAS_SLOTS * *n]  = 0;
	    bas[NPRIM_OF + BAS_SLOTS * *n]  = nexp;
	    bas[NCTR_OF  + BAS_SLOTS * *n]  = 1;
	    bas[PTR_EXP  + BAS_SLOTS * *n]  = *off;
	    bas[PTR_COEFF+ BAS_SLOTS * *n]  = *off+nexp;
	  }
	  (*n)++;
	  if(bas!=NULL && env!=NULL){
	    bas[ATOM_OF  + BAS_SLOTS * *n]  = atom_of;
	    bas[ANG_OF   + BAS_SLOTS * *n]  = 1;
	    bas[NPRIM_OF + BAS_SLOTS * *n]  = nexp;
	    bas[NCTR_OF  + BAS_SLOTS * *n]  = 1;
	    bas[PTR_EXP  + BAS_SLOTS * *n]  = *off;
	    bas[PTR_COEFF+ BAS_SLOTS * *n]  = *off+2*nexp;
	  }
	  (*n)++;
	  for(int i=0;i<nexp;i++){
	    fin.getline(buf,MAXLINELENGTH);
	    num=sscanf(buf,"%lf %lf %lf",&exp,&coef,&pcoef);
	    if(num!=3) throw std::invalid_argument( strcat(buf,"exp scoef pcoef") );
	    if(bas!=NULL && env!=NULL){
	      env[*off ] = exp;
	      env[*off+nexp]=coef;
	      env[*off+2*nexp]=pcoef;
	    }
	    (*off)++;
	  }	  
	  (*off)+=2*nexp;
	}	
      }else{
	if(found)break;
	for(int i=0;i<nexp;i++){
	  fin.getline(buf,MAXLINELENGTH);
	}
      }
    }
  }
  fin.close();
  return(*n);
}


// converts "SPDFGHI" into 01234567
int symlabel2ang(char s){
  char slabel[]="SPDFGHI";
  for(int i=0;i<7;i++){
    if(slabel[i]==s)return(i);
  }
}



ostream& operator<< (ostream &out, const basis_set & b){
  int i,j,n;
  out << "#natm=" << b.natm << endl;
  for(i=0;i<b.natm;i++){
    out << "#atm=" << i <<" "<< "Z=" << b.atm[CHARGE_OF + ATM_SLOTS * i] << endl;
    out << "#x=" << b.env[b.atm[PTR_COORD + ATM_SLOTS * i]+0] << endl;
    out << "#y=" << b.env[b.atm[PTR_COORD + ATM_SLOTS * i]+1] << endl;
    out << "#z=" << b.env[b.atm[PTR_COORD + ATM_SLOTS * i]+2] << endl;    
  }
  out << "#nbas=" << b.nbas << endl;
  for(n=0;n<b.nbas;n++){
    out << "#shell " << n << endl;
    out << "#on atom " << b.bas[ATOM_OF  + BAS_SLOTS * n] << endl;
    out << "#ang "     << b.bas[ANG_OF   + BAS_SLOTS * n] << endl;
    for(j=0;j<b.bas[NPRIM_OF + BAS_SLOTS * n];j++){
      out << "# "
	   << b.env[b.bas[PTR_EXP  + BAS_SLOTS * n]+j] << " "
	   << b.env[b.bas[PTR_COEFF+ BAS_SLOTS * n]+j] << endl; 
    }    
  }
  return(out);
}
