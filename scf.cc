#include <cstring>
#include <iostream>

#include <stdio.h>
#include <cmath>
#include "molecule.h"
#include "basis_set.h"
#include "matrix.h"
#include "scf.h"
using namespace std;



scf::scf(const molecule &arg_mol, const basis_set &arg_b) :
		mol(arg_mol), b(arg_b) {
	cout << "Entering SCF" << endl;
	h0 = new double[b.nbf * b.nbf];
	f = new double[b.nbf * b.nbf];
	c = new double[b.nbf * b.nbf];
	dmat = new double[b.nbf * b.nbf];
	oe = new double[b.nbf];

	double * temp = new double[b.nbf * b.nbf];
	double * fold=new double[b.nbf * b.nbf];

	for (int i = 0; i < b.nbf * b.nbf; i++){
		h0[i] = b.nuc[i] + b.kin[i];
		dmat[i]=0.;
		c[i]=0.;
	}
	double old_energy=0.;
	for(int it=0;it<40;it++){

	  construct_fmat(f, h0, dmat, b.teint);
	  print_iteration();
	  for (int i = 0; i < b.nbf * b.nbf; i++){
	    temp[i]=f[i];
	  }
	  transform(temp, b.olap_inv_sqrt, b.nbf);
	  eigval(temp, oe, b.nbf);
	  mmult(b.olap_inv_sqrt, temp, c, b.nbf);
	  
	  construct_dmat(dmat, c, mol.occ);
	  
	  energy=total_energy(dmat,h0, f);
	  if((energy-old_energy)*(energy-old_energy) < 1e-8 && it >2)break;
	}

	
	
	delete[] fold;
	delete[] temp;
}

scf::~scf() {
	delete[] c;
	delete[] dmat;

	delete[] oe;
	delete[] h0;
	delete[] f;
}


double scf::total_energy(double * dmat,double * h0, double * f) const{
  double * h0d=new double[b.nbf*b.nbf];
  double * fd=new double[b.nbf*b.nbf];
  mmult(f, dmat, fd, b.nbf);
  mmult(h0, dmat, h0d, b.nbf);
  double E=0.;
  for(int i=0;i<b.nbf;i++){
    E+=fd[i*b.nbf+i]+h0d[i*b.nbf+i];
  }
  delete[] h0d;
  delete[] fd;
  return(E);
}

int scf::print_iteration() {
	double * temp=new double[b.nbf*b.nbf];
	double * temp2=new double[b.nbf*b.nbf];

	printf("# Energy %e\n",energy);

	printf("# %10s","MO index:");
	for (int j = 0; j < b.nbf; j++) {
		printf(" %10d ", j);
	}
	printf("\n");
	printf("# %10s","orb. E:");
	for (int j = 0; j < b.nbf; j++) {
		printf(" %10.5f ", oe[j]);
	}
	printf("\n");
	//	printf("%10s","norm:");
	//	for (int j = 0; j < b.nbf; j++) {
	//		printf(" %10.7f ", norm(c + j * b.nbf, b.olap, b.nbf));
	//	}
	//	printf("\n");

	/*
	for (int i = 0; i < b.nbf; i++) {
		printf("%4s%5d:","bf",i);

		for (int j = 0; j < b.nbf; j++) {
			printf(" %10.5f ", c[j * b.nbf + i]);
		}
		printf("\n");
	}
	printf("f\n");
	for (int i = 0; i < b.nbf; i++) {
		printf("%4s%5d:","bf",i);
		for (int j = 0; j < b.nbf; j++) {
			printf(" %10.5f ", f[j * b.nbf + i]);
		}
		printf("\n");
	}
	
	
	printf("h0\n");
	for (int i = 0; i < b.nbf; i++) {
		printf("%4s%5d:","bf",i);
		for (int j = 0; j < b.nbf; j++) {
			printf(" %10.5f ", h0[j * b.nbf + i]);
		}
		printf("\n");
	}

	printf("c\n");
	for (int i = 0; i < b.nbf; i++) {
		printf("%4s%5d:","bf",i);
		for (int j = 0; j < b.nbf; j++) {
			printf(" %10.5f ", c[j * b.nbf + i]);
		}
		printf("\n");
	}

	printf("Dmat\n");
	for (int i = 0; i < b.nbf; i++) {
		printf("%4s%5d:","bf",i);
		for (int j = 0; j < b.nbf; j++) {
			printf(" %10.5f ", dmat[j * b.nbf + i]);
		}
		printf("\n");
	}
	*/
	/*	mmult(dmat, b.olap, temp, b.nbf);
	double trace=0.;
	for (int i = 0; i < b.nbf; i++) {
		trace+=temp[i*b.nbf+i];
	}
	printf("Tr Dmat*S=%f\n",trace);

	for (int i = 0; i < b.nbf*b.nbf; i++) temp2[i]=0;

	mmult(temp, dmat, temp2, b.nbf);
	printf("0.5*Dmat*S*0.5*Dmat-0.5 Dmat\n");
	for (int i = 0; i < b.nbf; i++) {
		printf("%4s%5d:","bf",i);
		for (int j = 0; j < b.nbf; j++) {
			printf(" %10.5f ", 0.25*temp2[j * b.nbf + i]-0.5*dmat[j * b.nbf + i]);
		}
		printf("\n");
	}

	*/
	/*
	for (int i = 0; i < b.nbf*b.nbf; i++) temp[i]=f[i];
	transform(temp, c, b.nbf);
	printf("Fmat\n");
	for (int i = 0; i < b.nbf; i++) {
		for (int j = 0; j < b.nbf; j++) {
			printf(" %10.5f ", temp[j * b.nbf + i]);
		}
		printf("\n");
	}


	for (int i = 0; i < b.nbf*b.nbf; i++) temp[i]=h0[i];
	transform(temp, c, b.nbf);
	printf("H0mat\n");
	for (int i = 0; i < b.nbf; i++) {
		for (int j = 0; j < b.nbf; j++) {
			printf(" %10.5f ", temp[j * b.nbf + i]);
		}
		printf("\n");
	}


	
	*/


	delete[] temp;

}

int scf::construct_dmat(double * dmat, double * c, const unsigned int * occ) {
	for (int mu = 0; mu < b.nbf; mu++) {
		for (int nu = 0; nu < b.nbf; nu++) {
			dmat[mu * b.nbf + nu] = 0.;
			for (int k = 0; k < b.nbf; k++) {
				dmat[mu * b.nbf + nu] += 0.5*occ[k] * c[k * b.nbf + mu] * c[k * b.nbf + nu];
			}
			dmat[nu * b.nbf + mu]=dmat[mu * b.nbf + nu];

		}
	}
}

void contribution_4ind(int i, int n, double t, double* f, double* dmat) {
	// 4 indices equal i,j,k,l
	f[i * n + i] += dmat[i * n + i] * t * 0.5;
}

void contribution_3ind(int i, int l, int n, double t, double* f, double* dmat) {
	// 3 indices equal i,j,k
	f[i * n + i] += dmat[i * n + l] * t;
	f[i * n + i] += dmat[l * n + i] * t;
	f[i * n + l] += dmat[i * n + i] * t;
	f[l * n + i] += dmat[i * n + i] * t;
	f[i * n + i] += -dmat[i * n + l] * t * 0.5;
	f[i * n + i] += -dmat[l * n + i] * t * 0.5;
	f[i * n + l] += -dmat[i * n + i] * t * 0.5;
	f[l * n + i] += -dmat[i * n + i] * t * 0.5;
}

void contribution_2indpairA(int i, int k, int n, double t, double* f,
		double* dmat) {
	// 2 indices equal i=j k=l
	f[i * n + i] += dmat[k * n + k] * t;
	f[k * n + k] += dmat[i * n + i] * t;
	f[i * n + k] += -dmat[i * n + k] * t * 0.5;
	f[k * n + i] += -dmat[k * n + i] * t * 0.5;
}

void contribution_2indpairB(int i, int j, int n, double t, double* f,
		double* dmat) {
	// 2 indices equal i=k j=l
	f[i * n + j] += dmat[i * n + j] * t;
	f[i * n + j] += dmat[j * n + i] * t;
	f[j * n + i] += dmat[i * n + j] * t;
	f[j * n + i] += dmat[j * n + i] * t;
	f[i * n + i] += -dmat[j * n + j] * t * 0.5;
	f[i * n + j] += -dmat[j * n + i] * t * 0.5;
	f[j * n + i] += -dmat[i * n + j] * t * 0.5;
	f[j * n + j] += -dmat[i * n + i] * t * 0.5;
}

void contribution_2indA(int i, int k, int l, int n, double t, double* f,
		double* dmat) {
//	return;
	// 2 indices equal i=j
	f[i * n + i] += dmat[k * n + l] * t;
	f[i * n + i] += dmat[l * n + k] * t;
	f[k * n + l] += dmat[i * n + i] * t;
	f[l * n + k] += dmat[i * n + i] * t;
	f[i * n + k] += -dmat[i * n + l] * t * 0.5;
	f[i * n + l] += -dmat[i * n + k] * t * 0.5;
	f[k * n + i] += -dmat[l * n + i] * t * 0.5;
	f[l * n + i] += -dmat[k * n + i] * t * 0.5;
}

void contribution_2indB(int i, int j, int l, int n, double t, double* f,
		double* dmat) {
//	return;
	// 2 indices equal i=k (j!=l)
	f[i * n + j] += dmat[i * n + l] * t;
	f[j * n + i] += dmat[i * n + l] * t;
	f[i * n + j] += dmat[l * n + i] * t;
	f[j * n + i] += dmat[l * n + i] * t;
	f[i * n + l] += dmat[i * n + j] * t;
	f[i * n + l] += dmat[j * n + i] * t;
	f[l * n + i] += dmat[i * n + j] * t;
	f[l * n + i] += dmat[j * n + i] * t;

	f[i * n + i] += -dmat[i * n + l] * t * 0.5;
	f[j * n + i] += -dmat[i * n + l] * t * 0.5;
	f[i * n + l] += -dmat[j * n + i] * t * 0.5;
	f[j * n + l] += -dmat[i * n + i] * t * 0.5;
	f[i * n + i] += -dmat[l * n + j] * t * 0.5;
	f[i * n + j] += -dmat[l * n + i] * t * 0.5;
	f[l * n + i] += -dmat[i * n + j] * t * 0.5;
	f[l * n + j] += -dmat[i * n + i] * t * 0.5;
}

void contribution_noind(int i, int j, int k, int l, int n, double t, double* f,
		double* dmat) {
//	return;
	// 2 indices equal i!=k (j!=l)
	f[i * n + j] += dmat[k * n + l] * t;
	f[j * n + i] += dmat[k * n + l] * t;
	f[i * n + j] += dmat[l * n + k] * t;
	f[j * n + i] += dmat[l * n + k] * t;
	f[k * n + l] += dmat[i * n + j] * t;
	f[k * n + l] += dmat[j * n + i] * t;
	f[l * n + k] += dmat[i * n + j] * t;
	f[l * n + k] += dmat[j * n + i] * t;

	f[i * n + k] += -dmat[j * n + l] * t * 0.5;
	f[j * n + k] += -dmat[i * n + l] * t * 0.5;
	f[i * n + l] += -dmat[j * n + k] * t * 0.5;
	f[j * n + l] += -dmat[i * n + k] * t * 0.5;
	f[k * n + i] += -dmat[l * n + j] * t * 0.5;
	f[k * n + j] += -dmat[l * n + i] * t * 0.5;
	f[l * n + i] += -dmat[k * n + j] * t * 0.5;
	f[l * n + j] += -dmat[k * n + i] * t * 0.5;
}

void scf::construct_fmat(double * f,double * h0, double * dmat, double * teint) {
  int n=b.nbf;
  for (int i = 0; i < b.nbf * b.nbf; i++)
    f[i] = h0[i];
  int ij, kl;
  int ijkl;
  int il, jk;
  int iljk;
  /*
    for(int i=0;i<n;i++) {
    for (int j = 0; j <= i; j++) {
    ij = j + i * (i - 1) / 2;
    for (int k = 0; k < n; k++) {
    for (int l = 0; l <= k; l++) {
    kl = l + k * (k - 1) / 2;
    if (ij < kl)
    continue;
    ijkl = kl + ij * (ij - 1) / 2;
    double t = teint[ijkl];
    
    if (i == j && k == l && i == k) {
    // 4 indices equal
    contribution_4ind(i, n, t, f, dmat);
    } else if (i == j && j == k) {
    // 3 indices equal i,j,k
    contribution_3ind(i, l, n, t, f, dmat);
    } else if (i == k && j == l) {
    // 3 indices equal i,j,l
    contribution_3ind(i, k, n, t, f, dmat);
    } else if (i == k && k == l) {
    // 3 indices equal i,k,l
    contribution_3ind(i, j, n, t, f, dmat);
    } else if (i == j && k == l) {
    // two indices pairwise equal ij and kl
    contribution_2indpairA(i, k, n, t, f, dmat);
    } else if (i == k && j == l) {
    // two indices pairwise equal ik and jl
    contribution_2indpairB(i, j, n, t, f, dmat);
    } else if (i == j) {
    // two indices  equal i=j
    contribution_2indA(i, k, l, n, t, f, dmat);
    } else if (k == l) {
    // two indices  equal k=l
    contribution_2indA(k, i, j, n, t, f, dmat);
    } else if (i == k) {
    // two indices  equal i=k
    contribution_2indB(i, j, l, n, t, f, dmat);
    } else if (j == l) {
    // two indices  equal j=l
    contribution_2indB(j, i, k, n, t, f, dmat);
    } else {
    contribution_noind(i, j, k, l, n, t, f, dmat);
    }//if
    }//l
    }//k
    }//j
    }//i
  */
  for (int i = 0; i < n; i++){
    for (int j = 0; j < n ; j++){
      if(i<j){ij = j + i * (i + 1)/ 2;}
      else{ij = i + j * (j + 1)/2;}
      ij=b.packed_index(i,j);
      
      for (int k = 0; k < n ; k++){
	for (int l = 0; l < n ; l++){
	  if(k>l){kl = l + k * (k + 1) / 2;}
	  else{kl = l + l * (l + 1) / 2;}
	  kl=b.packed_index(k,l);
	  if(i<l){il = l + i * (i + 1) / 2;}
	  else{il = i + l * (l + 1) / 2;}
	  il=b.packed_index(i,l);
	  if(k<j){jk = j + k * (k + 1) / 2;}
	  else{jk = k + j * (j + 1) / 2;}
	  jk=b.packed_index(j,k);
	  
	  if(kl<ij){ijkl = ij + kl * (kl + 1) / 2;}
	  else{ijkl = kl + ij * (ij + 1) / 2;}
	  if(il<jk){iljk = jk + il * (il + 1) / 2;}
	  else{iljk = il + jk * (jk + 1) / 2;}
	  ijkl=b.packed_index(ij,kl);
	  iljk=b.packed_index(il,jk);
	  
	  f[i*n+j]+=dmat[k*n+l]*(2.*teint[ijkl]-teint[iljk]);
	  // if(fabs(dmat[k*n+l])<1e-3)continue;
	  // printf("F[%d,%d]+= D(%d,%d) [ (%d %d|%d %d) - 0.5 (%d %d|%d %d) ]= %f [%f -0.5 %f]\n",
	  // 	 i,j,
	  // 	 k,l,
	  // 	 i,j,k,l,
	  // 	 i,l,j,k,
	  // 	 dmat[k*n+l],
	  // 	 teint[ijkl],
	  // 	 teint[iljk]
	  // 	 );
	  
	}
      }
    }
  }
	
}
