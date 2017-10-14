#include <cstring>
#include <iostream>

#include <stdio.h>

extern "C" {
int dsyev_(char *jobz, char *uplo, int *n, double *a, int *lda, double *w,
		double *work, int *lwork, int *info);
}

#include "molecule.h"
#include "basis_set.h"
#include "matrix.h"
#include "scf.h"
using namespace std;

scf::scf(const molecule &arg_mol, const basis_set &arg_b): mol(arg_mol) ,b(arg_b) {
	cout << "Entering SCF" << endl;
	h0 = new double[b.nbf * b.nbf];
	f = new double[b.nbf * b.nbf];
	c = new double[b.nbf * b.nbf];
	oe = new double[b.nbf];

	for (int i = 0; i < b.nbf * b.nbf; i++)
		h0[i] = b.nuc[i] + b.kin[i];
	memcpy(c,h0,sizeof(double)*b.nbf*b.nbf);
	eigval(c, oe, b.nbf);
	print_iteration();
}

scf::~scf() {
	delete[] c;
	delete[] oe;
	delete[] h0;
	delete[] f;
}


int scf::print_iteration() {
	printf("Coefficients:\n");
	for(int i=0;i<b.nbf;i++){
		for(int j=0;j<b.nbf;j++){
			printf(" %f ",c[i*b.nbf+j]);
		}
		printf("\n");
	}
}
