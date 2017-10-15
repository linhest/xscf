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



void mmult(double *a,double *b,double *c,int n){
/*	int i,j,k;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			c[i*n+j]=0.;
			for(k=0;k<n;k++){
				c[i*n+j]+=a[i*n+k]*b[k*n+j];
			}
		}
	}*/
	char transa='N';
	char transb='N';
	double alpha=1, beta=0.;
	dgemm_(&transa,&transb,&n,&n,&n,&alpha,a,&n,b,&n,&beta,c,&n);
}


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
		printf("%s\n",text);
	}else if(info>0){
		char text[100];
		sprintf(text,"convergence %d",info);
		throw std::invalid_argument( text );
	}
	return((int)info);
}
