
#include <deque>
#include <cstdio>
#include <cmath>
#include "diis.h"
#include "matrix.h"

using namespace std;

diis_set::diis_set(const int nsave, const int dim): nsave(nsave), dim(dim) {
  fSave.clear();
  errSave.clear();
  printf("# DIIS init\n");
}




double * diis_set::diis_interpolate (double* err, double* fmat)
{
  double * err_tmp=new double[dim];
  double * fmat_tmp=new double[dim];
  for(int i=0;i<dim;i++)err_tmp[i]=err[i];
  for(int i=0;i<dim;i++)fmat_tmp[i]=fmat[i];

  if(fSave.size()>=nsave){
      // remove oldest and add error matrix / F matrix to list
      delete[] fSave[0];
      delete[] errSave[0];
      fSave.pop_front();
      errSave.pop_front();
  }
  fSave.push_back(fmat_tmp);
  errSave.push_back(err_tmp);

  int N=(int)fSave.size();

  double * bmatrix=new double[(N+1)*(N+1)];
  double * rhs=new double[(N+1)];
  for(int i=0;i<N;i++){
      for(int j=i;j<N;j++){
	  double t=0;
//	  for(int k=0;k<sqrt(dim);k++){
//	      for(int l=0;l<sqrt(dim);l++)
//		t+=errSave[i][k*sqrt_dim+l]*errSave[j][l*sqrt_dim+k];
//	  }
	  for(int k=0;k<dim;k++)
	    t+=errSave[i][k]*errSave[j][k];
	  bmatrix[i*(N+1)+j]=t;
	  bmatrix[j*(N+1)+i]=t;
      }
  }
  for(int i=0;i<N;i++){
    bmatrix[i*(N+1)+N]=-1.;
    bmatrix[N*(N+1)+i]=-1.;
    rhs[i]=0.;
  }
  bmatrix[N*(N+1)+N]=0.;
  rhs[N]=-1.;

  // rescale Bmatrix
  if(bmatrix[0]>1.){
      double factor=1./bmatrix[0];
      for(int i=0;i<N;i++){
	 for(int j=0;j<N;j++){
	     bmatrix[i*(N+1)+j]*=factor;
	 }
      }
  }

//  print_matrix(bmatrix,N+1,(char *)"B");

  double * coef = new double[N+1];

  int error;
  error=solve_AxB(bmatrix, rhs, coef, N+1);
  if(error==0){
    for(int k=0;k<dim;k++)
      fmat[k]=0;

//    printf("%d %f \n",N,coef[N]);
    if(N==1)coef[0]=1.;
    for(int i=0;i<N;i++){
//	printf("%d %f \n",i,coef[i]);
	for(int k=0;k<dim;k++)
	  fmat[k]+=coef[i]*fSave[i][k];
    }
  }
  else{
      for(int k=0;k<dim;k++)
	fmat[k]=fmat_tmp[k];

  }

  delete[] bmatrix;
  delete[] rhs;
  delete[] coef;
  return(NULL);
}


diis_set::~diis_set(){
  for(int i=0;i<fSave.size();i++){
      delete[] fSave[i];
      delete[] errSave[i];
  }
  fSave.clear();
  errSave.clear();
}
