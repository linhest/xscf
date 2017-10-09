#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

#include <libiwl/iwl.h>
#include <libipv1/ip_lib.h>

#include <libpsio/psio.h>
#include <libchkpt/chkpt.h>
#include <libciomr/libciomr.h>
#include <psifiles.h>

#include "types.h"

#define MAXFACT 100

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))

extern "C"{
    FILE *infile, *outfile;
    char *psi_file_prefix;
}

extern "C"{ 
const char *gprgid() { const char *prgid = "MY"; return(prgid); }
}

namespace psi{
  namespace my{

  /* Function prototypes */
  void init_io(int argc, char *argv[]);
  void exit_io(void);
  int cc2unit(char *);

  }
} /* Namespace psi */
using namespace psi::my;


extern "C" void read_psi3();
extern "C" void psi_done();
#define IOFF_MAX 50604
static void init_ioff(void);
static int * ioff;
static void transformteints(double ** C, double * tsoints, int nso,int nmo,double **result);
static void transformteints2(double ** C, double * tsoints, int nso,int nmo,double **result);
static void transformH0(double ** C, double * T,double *V, int nso,int nmo,double **result);

int nmo;
int nso;
int nao;
double * oints,*tints,*zints;
char ** felement;
double * evals;
double enuc;
int nelec;
int ms2;
int * docc;
int * socc;
int nirreps;
int nsym;
int * sopi;
int rottype;
char * sym_label;
char * wfn;
double ** C;
double ** S;
double ** Cao;
double ** uso2ao;

double * Sao;

int *stype, *snuc, *snumg, *sprim,*sloc;
int nshell;
int natom;
int nprim;
int maxam;
int * l_length;

int * irrepmo;

char * llabel="spdfgh";

extern "C" double get_twoel(int i, int j, int k, int l);
extern "C" double get_onel(int i, int j);


void read_psi3(){
  int unit=72;
  int i,j,k,l;

  init_io(0, NULL);
  init_ioff();
  chkpt_init(PSIO_OPEN_OLD);

  if((chkpt_rd_ref() != 0) && (chkpt_rd_ref() != 2)){
    printf("only RHF / ROHF supported!\n");
    exit(2);
  }

  char * label=chkpt_rd_label();
  printf("#label: %s\n",label);
  ip_string("WFN", &wfn, 0);

  
  nsym=chkpt_rd_nsymhf();
  nirreps=chkpt_rd_nirreps();
  rottype=chkpt_rd_rottype();
  // rottype==3 linear molecule
  sym_label=chkpt_rd_sym_label();

  sopi=chkpt_rd_sopi();
  nso=chkpt_rd_nso();
  nao=chkpt_rd_nao();
  nmo=chkpt_rd_nmo();
  irrepmo=(int *)malloc(sizeof(int)*(nmo));
  
  int o;
  for(o=0,i=0;i<nirreps;i++){
    int n;
    for(n=0;n<sopi[i];n++,o++){
      irrepmo[o]=i;
    }
  }

  
  docc=chkpt_rd_clsdpi();
  socc=chkpt_rd_openpi();
  nelec=0;
  ms2=0;
  for(i=0;i<nsym;i++){
    nelec+=socc[i]+2*docc[i];
    ms2+=socc[i]*2;
  }


  int Nalpha=docc[0]+socc[0];
  int Nbeta=docc[0];
  sloc=chkpt_rd_sloc();   

  
  C=chkpt_rd_scf();
  evals=chkpt_rd_evals();



  // if we use special options we reorder the orbitals back 
  int * moorder=NULL;
  char * reorder=NULL;
  if(ip_exist("MOORDER",0) && ip_exist("REORDER",0) ) {
    
    ip_string("REORDER",&reorder,0);
    if(strcmp(reorder,"AFTER")==0){
      int size;
      ip_count("MOORDER",&size,0);
      if(size!=nmo){
	printf("Length of MOORDER != nmo\n");
	exit(-1);
      }
      
      moorder=(int *)malloc(sizeof(int)*size);
      for(i=0; i < size ; i++) 
	ip_data("MOORDER","%d",&(moorder[i]),1,i);
      
      printf("orbitals have been reordered: ordering back!\n");
      double ** Cold=(double **)malloc(sizeof(double*)*nso);
      double * oldevals=(double*)malloc(sizeof(double)*nmo);
      int mu;
      for(mu=0;mu<nso;mu++){
	Cold[mu]=(double*)malloc(sizeof(double)*nmo);
	for(i=0;i<nmo;i++){
	  Cold[mu][i]=C[mu][i];
	}
      }
      for(i=0;i<nmo;i++)
	oldevals[i]=evals[i];
      for(i=0;i<nmo;i++){
	evals[i]=oldevals[moorder[i]];
	printf("MO %d to MO %d\n",i,moorder[i]);	        
      }
      for(mu=0;mu<nso;mu++){
	for(i=0;i<nmo;i++){
	  C[mu][i]=Cold[mu][moorder[i]];
	}
      }
      
      for(i=0;i<nmo;i++)
	free(Cold[i]);
      free(Cold);
      free(oldevals);

    }
  }


  uso2ao=chkpt_rd_usotao();
  Cao = block_matrix(nao,nmo);
  mmult(uso2ao,1,C,0,Cao,0,nao,nso,nmo,0);


  int * sloc=chkpt_rd_sloc_new();   
  maxam=chkpt_rd_max_am();

  l_length=(int*)malloc(sizeof(int )*(maxam+1));
  l_length[0]=1;
  for(l=0;l<=maxam;l++){    
    if(l) l_length[l] = l_length[l-1] + l + 1;
  }

  snuc=chkpt_rd_snuc();  
  stype=chkpt_rd_stype();
  snumg=chkpt_rd_snumg();
  sprim=chkpt_rd_sprim();
  nprim=chkpt_rd_nprim();


  nshell=chkpt_rd_nshell();
  natom=chkpt_rd_natom();
  double ** geom =chkpt_rd_geom();
  double * zvals=chkpt_rd_zvals();
  
  felement=chkpt_rd_felement();
  
  printf("#Atoms:\n");
  for(i=0;i <natom;i++){
    printf("#%d %s\n",i,felement[i]);
  }


  int naoints=nao*(nao+1)/2;
  Sao=(double *)malloc(sizeof(double)*naoints);  
  iwl_rdone(PSIF_OEI,PSIF_AO_S,Sao,naoints,0,0,stdout);


  for(i=0;i <nmo;i++){
    double max=0.5;
    int center=-1;
    int thisshell=-1;
    int shell=0;
    int nucleus=0;
    for(j=0;j<nso;j++){

      if(shell+1 < nshell && sloc[shell+1]-1==j)shell++;
      //      printf("\t BF %d belongs to shell %d\n",j,shell);
      
      if(fabs(C[j][i]) > max){
	max=fabs(C[j][i]);
	center=j;
	thisshell=shell;
	nucleus=snuc[shell]-1;
      }
    }
    if(center > -1){
      //      printf("MO: %d nucleus: %d  %f\n",i,nucleus,max);
    }else{
      //      printf("MO %d is delocalized\n",i);
    }
  }




  enuc=chkpt_rd_enuc();

  int n_one=((nmo+1)*nmo)/2;
  int n_two=((n_one+1)*n_one)/2;

  oints = init_array(n_one);

  tints=init_array(n_two);
  printf("# %d twoelectron integrals\n",n_two);

  // special case with option CASSCF we do the transformation 
  // on our own
  if( strcmp(wfn,"CASSCF")==0 || moorder!=NULL){
    int n_one_nso=((nso+1)*nso)/2;
    int n_two_nso=((n_one_nso+1)*n_one_nso)/2;

    printf("reading H0 (SOBASIS)\n");
    double * Tints = init_array(n_one_nso);
    double * Vints = init_array(n_one_nso);
    iwl_rdone(PSIF_OEI,PSIF_SO_T , Tints, n_one_nso,
	      0, 0, stdout);
    iwl_rdone(PSIF_OEI,PSIF_SO_V , Vints, n_one_nso,
	      0, 0, stdout);
    printf("transforming H0 to MOBASIS\n");
    transformH0(C, Tints,Vints, nso,nmo,&oints);
    free(Tints); free(Vints);


    printf("reading TEI (SOBASIS)\n");
    double * tsoints=NULL;
    tsoints=init_array(n_two_nso);
    iwl_rdtwo(PSIF_SO_TEI, tsoints, ioff, nso,
	      0, 0, 
	      0, stdout);    

    printf("transforming %d TEI\n",n_two_nso);
    //    transformteints(C, tsoints, nso,nmo,&tints);
    transformteints2(C, tsoints, nso,nmo,&tints);
    free(tsoints);
    
  }else{

    iwl_rdone(PSIF_OEI,PSIF_MO_OEI, oints, n_one,
	      0, 0, stdout);
    
    iwl_rdtwo(PSIF_MO_TEI, tints, ioff, nmo,
	      0, 0, 
	      0, stdout);

  }

  chkpt_close();
  exit_io();


}



double get_onel(int i, int j)
{
   int ij ;
   double value ;

   if (i > j) {
      ij = ioff[i] + j;
      value = oints[ij] ;
      return(value) ;
      }
   else {
      ij = ioff[j] + i ;
      value = oints[ij] ;
      return(value) ;
      }
   return(oints[ij]) ;
}

double get_twoel(int i, int j, int k, int l)
{
   int ij, kl, ijkl ;

   ij = ioff[MAX0(i,j)] ;
   ij += MIN0(i,j) ;
   kl = ioff[MAX0(k,l)] ;
   kl += MIN0(k,l) ;
   ijkl = ioff[MAX0(ij,kl)] ;
   ijkl += MIN0(ij,kl) ;


   return(tints[ijkl]) ;
}


namespace psi {
  namespace my {

  void init_io(int argc, char *argv[])
  {
   char *progid;
  
   progid = (char *) malloc(strlen(gprgid())+2);
   sprintf(progid, ":%s",gprgid());
  
   //   psi_start(&infile,&outfile,&psi_file_prefix,argc-1, argv+1, 0);
   psi_start(&infile,&outfile,&psi_file_prefix,0, NULL, 0);
   ip_cwk_add(progid);
   free(progid);
  
   psio_init(); psio_ipv1_config();
  }

    void exit_io(void)
    {
      psio_done();
      psi_stop(infile,outfile,psi_file_prefix);
    }
    
    
    
  
  } /* Namespace tocprint */
} /* Namespace psi */



static void init_ioff(void)
{
   int i;

   /* set offsets for ij-type canonical ordering */
   ioff = (int *) malloc (IOFF_MAX * sizeof(int)) ;
   ioff[0] = 0;
   for (i = 1; i < IOFF_MAX ; i++) {
      ioff[i] = ioff[i-1] + i;
      }
}



static void transformH0(double ** C, double * T,double *V, int nso,int nmo,double **result){  
  int i,j,mu,nu;
  
  for(i=0;i<nmo;i++){
    for(j=0;j<=i;j++){
      int ij=i*(i+1)/2+j;
      (*result)[ij]=0.;

      for(mu=0;mu<nso;mu++){
	for(nu=0;nu<nso;nu++){
	  int munu;
	  if(mu>=nu)
	    munu=mu*(mu+1)/2+nu;
	  else
	    munu=nu*(nu+1)/2+mu;

	  (*result)[ij]+=C[mu][i]*C[nu][j]*(T[munu]+V[munu]);

	  
	}
      }      
    }
  }
}


void transformteints(double ** C, double * tsoints, int nso,int nmo,double **result){  
  int i,j,mu,nu;
  int k,l,ka,la;

  double * temp2=(double *)calloc(sizeof(double),nmo*(nmo+1)/2*nso*(nso+1)/2);

  //first contraction
#pragma omp parallel for schedule(dynamic) private(i,j,k,l,mu,nu,ka,la) shared(temp2,nso,nmo,result,C)
  for(i=0;i<nmo;i++){
    for(j=0; j<=i && j <nmo;j++){
      int ij=i*(i+1)/2+j;
	  
      for(mu=0;mu<nso;mu++){
	for(nu=0;nu<=mu &&nu<nso;nu++){
	  int munu=mu*(mu+1)/2+nu;
	
	  temp2[ij*nso*(nso+1)/2+munu]=0.;	      

          // contraction
	  for(ka=0;ka<nso;ka++){
	    for(la=0;la<nso;la++){
	      int kala;
	      if(ka>la)
		kala=ka*(ka+1)/2+la;
	      else
		kala=la*(la+1)/2+ka;
	      
	      int munukala;
	      if(munu>kala)
		munukala=munu*(munu+1)/2+kala;
	      else
		    munukala=kala*(kala+1)/2+munu;		  
	      
	      temp2[ij*nso*(nso+1)/2+munu]+=C[ka][i]*C[la][j]*tsoints[munukala];
	      
	    }//la
	  }//ka   
	  
	}//nu
      }//mu
      
    }//j
  }//i


  // second contraction
#pragma omp parallel for schedule(dynamic) private(i,j,k,l,mu,nu,ka,la) shared(temp2,nso,nmo,result,C)
  for(i=0;i<nmo;i++){
    for(j=0; j<=i && j <nmo;j++){
      int ij=i*(i+1)/2+j;

      for(k=0;k<nmo;k++){
	for(l=0; l<=k && l <nmo;l++){
	  int kl=k*(k+1)/2+l;
	  
	  if(kl>ij)continue;
	  int ijkl=ij*(ij+1)/2+kl;


	  (*result)[ijkl]=0.;
	  
	  
	  // contraction
	  for(ka=0;ka<nso;ka++){
	    for(la=0;la<nso;la++){
	      int kala;
	      if(ka>la)
		kala=ka*(ka+1)/2+la;
	      else
		kala=la*(la+1)/2+ka;
	      
	      
	      
	      (*result)[ijkl]+=C[ka][k]*C[la][l]* temp2[ij*nso*(nso+1)/2+kala];
	      
	    }//la
	  }//ka   
	  
	  //	  printf("%d %d %d %d =%e\n",i+1,j+1,k+1,l+1,get_twoel(i, j,k,l));
	  
	}//l
      }//k
      
    }//j
  }//i
  
  free(temp2);

  

}





void transformteints2(double ** C, double * tsoints, int nso,int nmo,double **result){  
  int i,j,mu,nu;
  int k,l,ka,la;

  int nmo2=nmo*nmo;
  int nso2=nso*nmo;
  int nmonso=nmo*nso;
  int nmo2_h=nmo*(nmo+1)/2;
  int nso2_h=nso*(nso+1)/2;

  double * temp2=(double *)calloc(sizeof(double),nmo2_h*nso2_h);
  double * temp=(double *)calloc(sizeof(double),nmonso*nso2_h);

  //  double * temp3=(double *)calloc(sizeof(double),nmo*nmo*nmo*nmo);
  //double * temp4=(double *)calloc(sizeof(double),nmo*nmo*nmo*nmo);


  //first contraction  (mu nu ka l) = sum_la C_la,l (mu nu ka la)
#pragma omp parallel for schedule(dynamic) private(i,j,k,l,mu,nu,ka,la) shared(temp,temp2,nso,nmo,result,C)
  for(mu=0;mu<nso;mu++){
    for(nu=0; nu<=mu ;nu++){
      int munu=mu*(mu+1)/2+nu;
      for(ka=0;ka<nso;ka++){
        for(l=0;l<nmo;l++){
          //  temp3[mu *nmo*nmo*nmo + nu *nmo*nmo + ka *nmo +l]=0;
          //         temp3[nu *nmo*nmo*nmo + mu *nmo*nmo + ka *nmo +l]=0;
          temp[(munu*nmonso)+ka*nmo+l]=0.;	      
          for(la=0;la<nso;la++){
            int kala;
            if(ka>la)
              kala=ka*(ka+1)/2+la;
            else
              kala=la*(la+1)/2+ka;
            int munukala;
            if(munu>kala)
              munukala=munu*(munu+1)/2+kala;
            else
              munukala=kala*(kala+1)/2+munu;		  
            temp[munu *nmonso + ka*nmo+l]+=C[la][l]*tsoints[munukala];	      
            //            temp3[mu *nmo*nmo*nmo + nu *nmo*nmo + ka *nmo +l]+=C[la][l]*tsoints[munukala];	       
            // if(mu!=nu)
            //  temp3[nu *nmo*nmo*nmo + mu *nmo*nmo + ka *nmo +l]+=C[la][l]*tsoints[munukala];	       
          }//la
        }//l
      }//ka
    }//nu
  }//mu

  //second contraction (mu nu k l) = sum_ka C_ka,k (mu nu ka l)
#pragma omp parallel for schedule(dynamic) private(i,j,k,l,mu,nu,ka,la) shared(temp,temp2,nso,nmo,result,C)
  for(mu=0;mu<nso;mu++){
    for(nu=0; nu<=mu ;nu++){
      int munu=mu*(mu+1)/2+nu;
      for(k=0;k<nmo;k++){
        for(l=0;l<=k;l++){
          int kl=k*(k+1)/2+l;
          temp2[munu * nmo2_h + kl]=0.;
          //          temp4[mu *nmo*nmo*nmo + nu *nmo*nmo + k *nmo +l]=0;
          //temp4[nu *nmo*nmo*nmo + mu *nmo*nmo + k *nmo +l]=0;
          for(ka=0;ka<nso;ka++){        
            temp2[munu * nmo2_h + kl]+=C[ka][k] * temp[munu * nmonso + ka*nmo+l];
            //            temp4[mu *nmo*nmo*nmo + nu *nmo*nmo + k *nmo +l]+=C[ka][k]*temp3[mu *nmo*nmo*nmo + nu *nmo*nmo + ka *nmo +l];
            //if(mu!=nu)
            //  temp4[nu *nmo*nmo*nmo + mu *nmo*nmo + k *nmo +l]+=C[ka][k]*temp3[mu *nmo*nmo*nmo + nu *nmo*nmo + ka *nmo +l];
          }//la
        }//i
      }//ka
    }//nu
  }//mu
  
  free(temp);
  temp=(double *)calloc(sizeof(double),nmonso*nmo2_h);

  // third contraction (k l mu j) = sum_nu C_nu,j (mu nu k l)
#pragma omp parallel for schedule(dynamic) private(i,j,k,l,mu,nu,ka,la) shared(temp,temp2,nso,nmo,result,C)
  for(k=0;k<nmo;k++){
    for(l=0;l<=k;l++){
      int kl=k*(k+1)/2+l;
      for(mu=0;mu<nso;mu++){ 
        for(j=0;j<nmo;j++){
          temp[kl*nmonso + mu*nmo + j]=0.;	      
          //temp3[k *nmo*nmo*nmo + l *nmo*nmo + mu *nmo +j]=0.;    
          //temp3[l *nmo*nmo*nmo + k *nmo*nmo + mu *nmo +j]=0.;    
          for(nu=0; nu<nso ;nu++){
            int munu;
            if(mu>nu)
              munu=mu*(mu+1)/2+nu;
            else
              munu=nu*(nu+1)/2+mu;
            temp[kl*nmonso + mu*nmo + j]+=C[nu][j] * temp2[munu * nmo2_h + kl];
            //temp3[k *nmo*nmo*nmo + l *nmo*nmo + mu *nmo +j]+=C[nu][j]* temp4[mu *nmo*nmo*nmo + nu *nmo*nmo + k *nmo +l];	      
            //if(k!=l)
            //              temp3[l *nmo*nmo*nmo + k *nmo*nmo + mu *nmo +j]+=C[nu][j]* temp4[mu *nmo*nmo*nmo + nu *nmo*nmo + k *nmo +l];	      
          }// nu
        }// l 
      }// mu
    }// j
  }// i

  // last contraction (ijkl) = (k l i j) = sum_mu C_mu,i (k l mu j)
#pragma omp parallel for schedule(dynamic) private(i,j,k,l,mu,nu,ka,la) shared(temp,temp2,nso,nmo,result,C)
  for(k=0;k<nmo;k++){
    for(l=0;l<=k;l++){
      int kl=k*(k+1)/2+l;
      for(i=0;i<nmo;i++){ 
        for(j=0;j<=i;j++){
          int ij=i*(i+1)/2+j;
          if(kl>ij)continue;
	  int ijkl=ij*(ij+1)/2+kl;
          (*result)[ijkl]=0.;
          for(mu=0;mu<nso;mu++){
            (*result)[ijkl]+=C[mu][i] * temp[kl*nmonso + mu*nmo + j];
            //(*result)[ijkl]+=C[mu][i] *temp3[k *nmo*nmo*nmo + l *nmo*nmo + mu *nmo +j];
          }// mu
        }// j
      }// i
    }// l
  }// k


  free(temp);
  free(temp2);

}
