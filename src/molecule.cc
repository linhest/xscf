#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>

#include <vector>

#include "molecule.h"
#include "atom_labels.h"
using namespace std;

extern const char  atom_labels[][4];

#define MAXLINELENGTH 2000
molecule::molecule(char filename[]){
  ifstream fin;
  fin.open(filename);
  if(!fin.good()) {
	  throw std::invalid_argument( strcat(filename," not found ") );
  }
  while(!fin.eof()){
    char buf[MAXLINELENGTH];
    fin.getline(buf,MAXLINELENGTH);

    int l=strlen(buf);
    if(l==MAXLINELENGTH){
      throw std::invalid_argument( "line too long" );
    }
    if(strchr(buf,'#')!=NULL){
      char * a=strchr(buf,'#');
      *a='\0';
    }
    // concatenate lines if end with \\n
    while(l>2&&buf[l-1]=='\\'){
      fin.getline(buf+(l-1),MAXLINELENGTH-l);
      l=strlen(buf);
      if(l==MAXLINELENGTH)
	throw std::invalid_argument( "line too long" );      
    }

    
    char * token=strtok(buf,"=");
    if(token==NULL)continue;
    char * token2=strtok(NULL,"=");
    if(token2==NULL)continue;
    if(strcasecmp(token,"nref")==0){
      nref=parse_int(token2);
    }else if(strcasecmp(token,"basis_set")==0){
      strncpy(basis_set,token2,LEN_BASIS_SET);
    }else if(strcasecmp(token,"exlevels")==0){
      exlevels=parse_int(token2);
    }else if(strcasecmp(token,"occ")==0){
      occ=new unsigned int[nmo];
      parse_int_array_fixed_len(occ,nmo,token2,',');
    }else if(strcasecmp(token,"charge")==0){
      charge=parse_int(token2);
    }else if(strcasecmp(token,"sym")==0){
      target_sym=parse_int(token2);
    }else if(strcasecmp(token,"nmo")==0){
      nmo=parse_int(token2);
    }else if(strcasecmp(token,"ref")==0){
      vector<char *>* list=parse_string_array( token2, ';');
      nref=list->size();
      refs=new unsigned int[nref*nmo];
      int i;
      for(i=0;i<nref;i++){
	parse_int_array_fixed_len(&refs[i*nmo],nmo,(*list)[i],',');
      }
      delete list;
    }else if(strcasecmp(token,"geom")==0){
      char * ptr_start=token2;
      l=strlen(token)+strlen(token2);
      ptr_start=strchr(token2,'(')-1;
      while(ptr_start==NULL){
	fin.getline(buf+l-1,MAXLINELENGTH);	
	ptr_start=strchr(ptr_start,'(');
	l=strlen(token)+strlen(token2);
      }
      // concatenate lines until ')'
      while(strchr(ptr_start,')')==NULL){
	*(ptr_start+strlen(ptr_start))='\n';
	fin.getline(ptr_start+strlen(ptr_start)+1,MAXLINELENGTH-strlen(token)-strlen(token2));
       	l=strlen(token)+strlen(token2)+3;
	if(l>=MAXLINELENGTH)
	  throw std::invalid_argument( "line too long" );      
      }
      char * ptr_end=strchr(ptr_start,')');
      *(ptr_end+1)='\0';
      *(ptr_end)='\0';
      ptr_start=ptr_start+2;
      l=strlen(buf);
      if(l==MAXLINELENGTH)
	throw std::invalid_argument( "line too long" );
      vector< double* > G_list;
      vector<char*> label_list;
      int i;
      char *ptr, *ptr_next;
      for(i=0,ptr=ptr_start+1;ptr!=NULL;){
	char label[4];
	double g[3];
	int num;
	ptr_next=strchr(ptr,'\n');
	if(ptr_next!=NULL)
	  *(ptr_next)='\0';
	if(strlen(ptr)!=0){
	  num=sscanf(ptr,"%s %lf %lf %lf ",label,&g[0],&g[1],&g[2]);
	  if(num==4){
	    double * g_new = new double[3];
	    char * l_new = new char[4];
	    g_new[0]=g[0];
	    g_new[1]=g[1];
	    g_new[2]=g[2];
	    strncpy(l_new,label,4);
	    G_list.push_back(g_new);
	    label_list.push_back(l_new);
	    i++;
	  }else if(num>0){
	    throw std::invalid_argument( "error reading geom" );
	  }
	}
	if(ptr_next==NULL)break;
	ptr=ptr_next+1;
      }
      natm=G_list.size();
      geom=new double[natm*3];
      Z=new int[natm];
      atm_label=new char *[natm];
      for(i=0;i<natm;i++){
	geom[i*3+0]=G_list[i][0];
	geom[i*3+1]=G_list[i][1];
	geom[i*3+2]=G_list[i][2];
	atm_label[i]=label_list[i];
	for(int j=0;j<N_ATOM_LABELS;j++){
	  if(strncasecmp(atm_label[i],atom_labels[j],4)==0) Z[i]=j;
	}
      }
      G_list.empty();
      label_list.empty();
      
    }else{      
      throw std::invalid_argument( strcat(token," not recognized") );      
    }// check token
  }//!feof
  fin.close();
    
}// initialize


molecule::~molecule(){
  delete[] refs;
}

ostream& operator<< (ostream &out, const molecule & mol){
  int i,j;
  out << "#nmo="<< mol.nmo << endl;
  out << "#ms2="<< mol.ms2 << endl;
  out << "#nref="<< mol.nref << endl;
  out << "#basis_set="<< mol.basis_set << endl;
  out << "#Geometry=" << endl;
  for(i=0;i<mol.natm;i++){
    printf("# %d %d %4s %f %f %f\n",
	   i,
	   mol.Z[i],
	   mol.atm_label[i],
	   mol.geom[i*3+0],
	   mol.geom[i*3+1],
	   mol.geom[i*3+2]);
  }
  out << "#References=" << endl;
  for(i=0;i<mol.nref;i++){
    out << "#"<< i << ":" ;
    for(j=0;j<mol.nmo;j++){
      out << mol.refs[i*mol.nmo+j];
      if(j<mol.nmo-1){
	out << ",";
      }
    }
    out << endl;
  }
  return(out);
}

char *  molecule::get_basis_set(){
  return(&basis_set[0]);
}
int  molecule::get_natm(){
  return(natm);
}

int molecule::parse_int(char * string){
  int n;
  stringstream ss;
  ss << string;
  ss>>n;
  return(n);
}

vector<char *>* molecule::parse_string_array(char * buf, const char delim=','){
  vector<char *>* s = new vector<char *>();
  char * token;
  int n;
  for(n=0,token=strtok(buf,&delim); token!=NULL;token=strtok(NULL,&delim)){
    s->push_back(token);
   }
  return(s);
}

unsigned int * molecule::parse_int_array_fixed_len(unsigned int * s,  int nmo, char * buf, const char delim=','){
  char * token;
  int i;
  for(i=0,token=strtok(buf,&delim); token!=NULL && i <nmo;token=strtok(NULL,&delim),i++){
    s[i]=parse_int(token);
  }
  if(i==nmo && token !=NULL)
    throw std::invalid_argument( "array too long" );      
  for(;i<nmo;i++){
    s[i]=0;
  }
  
  return(s);
}

int * molecule::parse_int_array_fixed_len(int * s,  int nmo, char * buf, const char delim=','){
  char * token;
  int i;
  for(i=0,token=strtok(buf,&delim); token!=NULL && i <nmo;token=strtok(NULL,&delim),i++){
    s[i]=parse_int(token);
  }
  if(i==nmo && token !=NULL)
    throw std::invalid_argument( "array too long" );      
  for(;i<nmo;i++){
    s[i]=0;
  }
  return(s);
}



				    
