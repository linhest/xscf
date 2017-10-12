#include <cstring>
#include <iostream>

#include <stdio.h>
extern "C" {
int dsyev_(char *jobz, char *uplo, int *n, double *a,
	 int *lda, double *w, double *work, int *lwork,
	int *info);

}
using namespace std;


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
		sprintf(text,"argument %ld",-info);
		throw std::invalid_argument( text );
		printf("%s\n",text);
	}else if(info>0){
		char text[100];
		sprintf(text,"convergence %ld",info);
		throw std::invalid_argument( text );
	}
	return((int)info);
}
