#include <cstring>
#include <iostream>
#include <fstream>
#include <deque>

#include <stdio.h>
#include <cmath>
#include "diis.h"
#include "matrix.h"
#include "molecule.h"

#include "basis_set.h"

#include "scf.h"
using namespace std;



scf::scf(const molecule &mol, const basis_set &b) :
  mol(mol),
  b(b),
  diis(mol.n_diis,b.nbf*b.nbf){
	cout << "# Entering SCF" << endl;

	h0 = new double[b.nbf * b.nbf];
	f = new double[b.nbf * b.nbf];
	c = new double[b.nbf * b.nbf];
	dmat = new double[b.nbf * b.nbf];
	dmat_old = new double[b.nbf * b.nbf];
	oe = new double[b.nbf];

	double * temp = new double[b.nbf * b.nbf];
	double * fold=new double[b.nbf * b.nbf];

	for (int i = 0; i < b.nbf * b.nbf; i++){
		h0[i] = b.nuc[i] + b.kin[i];
		dmat[i]=0.;
		c[i]=0.;
	}
	if(mol.print>2)
	  print_matrix(h0,b.nbf,(char *)"H0 Matrix");

	if(mol.filename_guess_read[0]!='\0')
	  read_guess(mol.filename_guess_read);

	energy=0.;
	for(it=0;it<40;it++){
	  old_energy=energy;

	  for (int i = 0; i < b.nbf * b.nbf; i++)dmat_old[i]=dmat[i];
	  construct_dmat(dmat, c, mol.occ);

	  if(mol.print>2)
	    print_matrix(dmat,b.nbf,(char *)"D Matrix");

	  construct_fmat(f, h0, dmat, b.teint);




	  // DIIS Interpolation
	  if(it>0  && mol.n_diis>0 ){
	    double * err = new double[b.nbf*b.nbf];
	    diis_error(dmat, f, err);
	    diis.diis_interpolate (err, f);
	    delete[] err;
	  }

	  for (int i = 0; i < b.nbf * b.nbf; i++)temp[i]=f[i];
	  transform(temp, b.olap_inv_sqrt, b.nbf);
	  eigval(temp, oe, b.nbf);
	  mmult(b.olap_inv_sqrt, temp, c, b.nbf);


	  energy=total_energy(dmat,h0, f)+mol.Enuc;
	  max_dmat_diff=0.;
	  for(int i=0;i<b.nbf*b.nbf;i++)
	    if(max_dmat_diff < fabs(dmat[i]-dmat_old[i]))
	      max_dmat_diff=fabs(dmat[i]-dmat_old[i]);

	  print_iteration();


	  if( it >1
	      && fabs(energy-old_energy) < mol.energy_conv
	      && max_dmat_diff < mol.density_conv)
	    break;

	}

	if( it >1
	    && fabs(energy-old_energy) < mol.energy_conv
	    && max_dmat_diff < mol.density_conv){
	    print_result();
	}else{
	    printf("CONVERGENCE FAILED\n");
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

double scf::diis_error(double * dmat, double *f, double * err){
  int n=b.nbf;

  double * fd = new double[n*n];
  double * fds = new double[n*n];
  mmult(f, dmat, fd, n);
  mmult(fd,b.olap, fds, n);


  double * sd = new double[n*n];
  double * sdf = new double[n*n];
  mmult(b.olap, dmat, sd, n);
  mmult(sd, f, sdf, n);

  double t=0.;
  for(int i=0;i<n*n;i++){
      if(err!=NULL)err[i]=sdf[i]-fds[i];
      if(fabs(sdf[i]-fds[i])>t)t=fabs(sdf[i]-fds[i]);
  }

  delete[] sd;
  delete[] sdf;
  delete[] fd;
  delete[] fds;
  return(t);
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
	if(it==0)
	  printf("# %5s %15s %15s %15s %15s\n",
		 "Iter",
		 "Energy",
		 "Energy diff",
		 "Density diff",
		 "max |FDS-SDF|");
	else
	  printf("# %5d %15.5e %15.5e %15.5e %15.5e\n",
		 it,
		 energy,
		 energy-old_energy,
		 max_dmat_diff,
		 diis_error(dmat, f, NULL));
	if(mol.print>1){
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
	}
}


int scf::print_result() {

	printf("# Results:\n# Energy=%f \n# Delta E=%e Delta D=%e DIIS Error=%e\n",energy,
	       energy-old_energy,
	       max_dmat_diff,
	       diis_error(dmat, f, NULL));
	printf("# %10s","MO index:");
	for (int j = 0; j < b.nbf; j++) {
	    printf(" %10d ", j);
	}
	printf("\n");
	printf("# %10s","orb. E:");
	for (int j = 0; j < b.nbf; j++) {
	    printf(" %10.4f ", oe[j]);
	}
	printf("\n");

	for (int i = 0; i < b.nbf; i++) {
		printf("# %4s%5d:","bf",i);

		for (int j = 0; j < b.nbf; j++) {
			printf(" %10.7f ", c[j * b.nbf + i]);
		}
		printf("\n");
	}
}

int scf::save_guess(const char * filename){
  printf("saving to %s\n",filename);
  fstream fout(filename,ios::out|ios::binary|ios::trunc);
  if (!fout.good())
    throw std::invalid_argument(" cannot write ");
  fout.write((char *)c,sizeof(double)*b.nbf*b.nbf);
  fout.close();
  return(0);
}

int scf::read_guess(const char * filename){
  printf("reading from %s\n",filename);
  fstream fin(filename,ios::in|ios::binary);
  if (!fin.good())
    throw std::invalid_argument(" cannot read ");
  fin.read((char *)c,sizeof(double)*b.nbf*b.nbf);
  fin.close();
  return(0);
}



int scf::construct_dmat(double * dmat, double * c, const unsigned int * occ) {
	for (int mu = 0; mu < b.nbf; mu++) {
		for (int nu = mu; nu < b.nbf; nu++) {
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
      ij=b.packed_index(i,j);
      for (int k = 0; k < n; k++) {
	for (int l = 0; l <= k; l++) {
	  kl=b.packed_index(k,l);
	  if (ij < kl)
	    continue;
	  ijkl=b.packed_index(ij,kl);
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

      ij=b.packed_index(i,j);
      
      for (int k = 0; k < n ; k++){
	for (int l = 0; l < n ; l++){
	  kl=b.packed_index(k,l);
	  il=b.packed_index(i,l);
	  jk=b.packed_index(j,k);
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
