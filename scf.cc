#include <cstring>
#include <iostream>

#include <stdio.h>


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
	dmat = new double[b.nbf * b.nbf];

	oe = new double[b.nbf];

	for (int i = 0; i < b.nbf * b.nbf; i++)
		h0[i] = b.nuc[i] + b.kin[i];

	for (int i = 0; i < b.nbf * b.nbf; i++)
		f[i] = h0[i];
	transform(f,b.olap_inv_sqrt,b.nbf);
	memcpy(c,f,sizeof(double)*b.nbf*b.nbf);
	eigval(c, oe, b.nbf);
	construct_dmat(dmat, c, mol.occ );

	for(int i=0;i<b.nbf;i++){
		printf("%d %f %f \n",i,oe[i],norm(c+i*b.nbf,b.olap,b.nbf));
	}

	print_iteration();


}

scf::~scf() {
	delete[] c;
	delete[] dmat;

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

int scf::construct_dmat(double * dmat, double * c, const unsigned int * occ ){
	for(int mu=0;mu<b.nbf;mu++){
		for(int nu=0;nu<b.nbf;nu++){
			dmat[mu*b.nbf+nu]=0.;
			for(int k=0;k<b.nbf;k++){
				dmat[mu*b.nbf+nu]=occ[k]*c[k*b.nbf+mu]*c[k*b.nbf+nu];
			}

		}
	}
}

