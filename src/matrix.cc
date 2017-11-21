#include <cstring>
#include <iostream>

#include <stdio.h>
extern "C" {
int dsyev_(char *jobz, char *uplo, int *n, double *a,
	 int *lda, double *w, double *work, int *lwork,
	int *info);

int dsymm_(char *side, char *uplo, int *m, int *n,
	double *alpha, double *a, int *lda, double *b,
	int *ldb, double *beta, double *c, int *ldc);
int dgemm_(char *transa, char *transb, int *m, int *
	n, int *k, double *alpha, double *a, int *lda,
	double *b, int *ldb, double *beta, double *c, int
	*ldc);
int dgemv_(char *trans, int *m, int *n, double *
	alpha, double *a, int *lda, double *x, int *incx,
	double *beta, double *y, int *incy);
int dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv,
                double* b, int* ldb, int * info);

double ddot_(int *n, double *dx, int *incx, double *dy, int *incy);
}
using namespace std;



// a is rewritten as b^t a b
void transform(double *a,double *b,int n){
	double * temp=new double[n*n];
	char transa='N';
	char transb='N';
	double alpha=1, beta=0.;
	dgemm_(&transa,&transb,&n,&n,&n,&alpha,a,&n,b,&n,&beta,temp,&n);//temp=a*b
	transa='T';
	transb='N';
	dgemm_(&transa,&transb,&n,&n,&n,&alpha,b,&n,temp,&n,&beta,a,&n);//a=b^t*temp
	delete [] temp;
}

// return C_mu S_munu C_nu
double norm(double *C,double *S,int n){
	double * temp=new double[n];
	char transa='N';
	double alpha=1,beta=0.;
	int inc=1;
	dgemv_(&transa,&n,&n,&alpha,S,&n,C,&inc,&beta,temp,&inc);//temp=S*C
	double norm_val=ddot_(&n, C, &inc, temp, &inc);
	delete [] temp;
	return(norm_val);
}



// a is rewritten as c=a b
void mmult(double *a,double *b,double * c,int n){
	char transa='N';
	char transb='N';
	double alpha=1, beta=0.;
	dgemm_(&transa,&transb,&n,&n,&n,&alpha,a,&n,b,&n,&beta,c,&n);//c=a*b
}

// calculate eigenvectors + eigenvalues of real symmetric matrix of dimension m
int eigval(double * a, double * eval, int m){
	int n,lwork,info;
	n=m;
	char jobz='V';
	char uplo='L';
	double * work;
	double opt_lwork=0.0;
	int r;
	lwork=-1;
	r=dsyev_(&jobz, &uplo, &n, a, &n, eval, (double *)&opt_lwork, &lwork, &info);
	lwork=(int)opt_lwork;
	work=new double[lwork];
	r=dsyev_(&jobz, &uplo, &n, a, &n, eval, work, &lwork, &info);
	delete[] work;
	if(info<0){
		char text[100];
		sprintf(text,"argument %d",-info);
		throw std::invalid_argument( text );
	}else if(info>0){
		char text[100];
		sprintf(text,"convergence %d",info);
		throw std::invalid_argument( text );
	}
	return((int)info);
}

int solve_AxB(double * a, double *b, double * d, int n){
  int * ipiv=new int[n];
  int nrhs=1;
  int info;
  double * c=new double[n*n];
  for(int i=0;i<n*n;i++)c[i]=a[i];
  for(int i=0;i<n;i++)d[i]=b[i];
  dgesv_( &n, &nrhs, c, &n, ipiv,d, &n, &info );
  delete[] ipiv;
  delete[] c;
  return(info);
}


double trace(double *m, int n){
  double t=0.;
  for(int i=0;i<n;i++)
      t+=m[i*n+i];
  return(t);
}


void print_matrix(double *m, int n, char * description){
  printf("Matrix %s:\n",description);
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      printf("%f ",m[j*n+i]);
    }
    printf("\n");
  }

}
