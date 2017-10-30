#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstring>
#include <cmath>

#include "matrix.h"
#include "molecule.h"
#include "basis_set.h"
using namespace std;

extern "C" {
#include <cint.h>
  int cint1e_ovlp_sph(double *buf, int *shls, int *atm, int natm, int *bas,
		      int nbas, double *env);
  int cint1e_kin_sph(double *buf, int *shls, int *atm, int natm, int *bas,
		     int nbas, double *env);
  int cint1e_nuc_sph(double *buf, int *shls, int *atm, int natm, int *bas,
		     int nbas, double *env);

  int dgemm_(char *transa, char *transb, int *m, int * n, int *k, double *alpha, double *a, int *lda,
	     double *b, int *ldb, double *beta, double *c, int *ldc);
  int dgemv_(char *trans, int *m, int *n, double *
	     alpha, double *a, int *lda, double *x, int *incx,
	     double *beta, double *y, int *incy);
}
#define ENV_SIZE 10000
#define MAXLINELENGTH 2000

basis_set::~basis_set() {
  delete[] atm;
  delete[] bas;
  delete[] env;
}
basis_set::basis_set(molecule & mol) {
  char * filename = mol.get_basis_set();
  int off = PTR_ENV_START; // = 20
  int n = 0;
  for (int i = 0; i < mol.natm; i++) {
    read_basis_set_atom(filename, mol.atm_label[i], NULL, NULL, i, &off,
			&n);
    cout << mol.atm_label[i] << " " << n << endl;
  }
  nbas = n;
  natm = mol.natm;
  atm = new int[mol.natm * ATM_SLOTS];
  bas = new int[nbas * BAS_SLOTS];
  env = new double[ENV_SIZE];
  off = PTR_ENV_START; // = 20
  n = 0;
  for (int i = 0; i < mol.natm; i++) {
    atm[CHARGE_OF + ATM_SLOTS * i] = mol.Z[i];
    atm[PTR_COORD + ATM_SLOTS * i] = off;
    env[off + 0] = mol.geom[i * 3 + 0];
    env[off + 1] = mol.geom[i * 3 + 1];
    env[off + 2] = mol.geom[i * 3 + 2];
    off += 3;
  }
  for (int i = 0; i < mol.natm; i++) {
    read_basis_set_atom(filename, mol.atm_label[i], env, bas, i, &off, &n);
  }
  nbf = 0;
  lmax = 0;
  for (int n = 0; n < nbas; n++) {
    if (lmax < bas[ANG_OF + BAS_SLOTS * n])
      lmax = bas[ANG_OF + BAS_SLOTS * n];
  }
  first_bf_of_shell = new int[nbas];
  shells_on_atom = new int[natm * lmax];
  for (int n = 0; n < nbas; n++) {
    int l = bas[ANG_OF + BAS_SLOTS * n];
    int a = bas[ATOM_OF + BAS_SLOTS * n];
    first_bf_of_shell[n] = nbf;
    nbf += 2 * l + 1;
  }

  olap = calc_olap();
  kin = calc_kin();
  nuc = calc_nuc();
  teint = calc_teint();
  olap_inv_sqrt = inv_sqrt(olap);
}

// reads in basis set data
// from file filename
// for atom atm_label with index atom_of
// env,bas,off,n are the libcint arrays
int basis_set::read_basis_set_atom(char * filename, char * atm_label,
				   double * env, int * bas, int atom_of, int *off, int *n) {
  ifstream fin;
  bool found = false;
  int nshell_atm = 0;
  fin.open(filename);
  if (!fin.good())
    throw std::invalid_argument(strcat(filename, " not found "));
  while (!fin.eof()) {
    char alabel[4];
    char slabel[2];
    int nexp, nshell;
    double exp, norm, coef, pcoef;
    int ai;
    int num;
    char buf[MAXLINELENGTH];
    fin.getline(buf, MAXLINELENGTH);
    if (strncmp("****", buf, 4) == 0) {
      // new entry starts here
      fin.getline(buf, MAXLINELENGTH);
      if (fin.eof())
	break;
      num = sscanf(buf, "%s %d", alabel, &ai);
      if (num != 2)
	throw std::invalid_argument(strcat(buf, "alabel int"));
    } else {
      num = sscanf(buf, "%s %d %lf", slabel, &nexp, &norm);
      if (num != 3)
	throw std::invalid_argument(strcat(buf, "sym nexp norm"));
      if (strcasecmp(alabel, atm_label) == 0) {
	found = true;
	if (strlen(slabel) == 1) {
	  nshell_atm++;
	  if (bas != NULL && env != NULL) {
	    bas[ATOM_OF + BAS_SLOTS * *n] = atom_of;
	    bas[ANG_OF + BAS_SLOTS * *n] = symlabel2ang(slabel[0]);
	    bas[NPRIM_OF + BAS_SLOTS * *n] = nexp;
	    bas[NCTR_OF + BAS_SLOTS * *n] = 1;
	    bas[PTR_EXP + BAS_SLOTS * *n] = *off;
	    bas[PTR_COEFF + BAS_SLOTS * *n] = *off + nexp;
	  }
	  for (int i = 0; i < nexp; i++) {
	    fin.getline(buf, MAXLINELENGTH);
	    num = sscanf(buf, "%lf %lf ", &exp, &coef);
	    if (num != 2)
	      throw std::invalid_argument(strcat(buf, "exp coef"));
	    if (bas != NULL && env != NULL) {
	      env[*off] = exp;
	      env[*off + nexp] = coef
		* CINTgto_norm(bas[ANG_OF + BAS_SLOTS * (*n)],exp);
					    
	    }
	    (*off)++;
	  }
	  (*off) += nexp;
	  (*n)++;
	} else if (strlen(slabel) == 2 && strcmp(slabel, "SP") == 0) {
	  nshell_atm += 2;
	  if (bas != NULL && env != NULL) {
	    bas[ATOM_OF + BAS_SLOTS * *n] = atom_of;
	    bas[ANG_OF + BAS_SLOTS * *n] = 0;
	    bas[NPRIM_OF + BAS_SLOTS * *n] = nexp;
	    bas[NCTR_OF + BAS_SLOTS * *n] = 1;
	    bas[PTR_EXP + BAS_SLOTS * *n] = *off;
	    bas[PTR_COEFF + BAS_SLOTS * *n] = *off + nexp;
	  }
	  (*n)++;
	  if (bas != NULL && env != NULL) {
	    bas[ATOM_OF + BAS_SLOTS * *n] = atom_of;
	    bas[ANG_OF + BAS_SLOTS * *n] = 1;
	    bas[NPRIM_OF + BAS_SLOTS * *n] = nexp;
	    bas[NCTR_OF + BAS_SLOTS * *n] = 1;
	    bas[PTR_EXP + BAS_SLOTS * *n] = *off;
	    bas[PTR_COEFF + BAS_SLOTS * *n] = *off + 2 * nexp;
	  }
	  (*n)++;
	  for (int i = 0; i < nexp; i++) {
	    fin.getline(buf, MAXLINELENGTH);
	    num = sscanf(buf, "%lf %lf %lf", &exp, &coef, &pcoef);
	    if (num != 3)
	      throw std::invalid_argument(
					  strcat(buf, "exp scoef pcoef"));
	    if (bas != NULL && env != NULL) {
	      env[*off] = exp;
	      env[*off + nexp] = coef;
	      env[*off + 2 * nexp] = pcoef;
	    }
	    (*off)++;
	  }
	  (*off) += 2 * nexp;
	}
      } else {
	if (found)
	  break;
	for (int i = 0; i < nexp; i++) {
	  fin.getline(buf, MAXLINELENGTH);
	}
      }
    }
  }
  fin.close();
  // Normalize the basis
  /*
  if (bas != NULL && env != NULL) {
    int sh;
    int shls[2];
    for (sh = 0; sh < (*n); sh++) {
      int nl = 2 * bas[ANG_OF + BAS_SLOTS * sh] + 1;
      double * buf = new double[nl * nl];
      shls[0] = sh;
      shls[1] = sh;
      cint1e_ovlp_sph(buf, shls, atm, natm, bas, nbas, env);
      for(int i=0;i< bas[NPRIM_OF + BAS_SLOTS * sh];i++){
	env[bas[PTR_COEFF + BAS_SLOTS * sh]+i ] /= sqrt(buf[0]);
      }
      delete[] buf;
    }
  }
  */
	
  return (*n);
}

// converts "SPDFGHI" into 01234567
int symlabel2ang(char s) {
  char slabel[] = "SPDFGHI";
  for (int i = 0; i < 7; i++) {
    if (slabel[i] == s)
      return (i);
  }
}

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))

int basis_set::packed_index (int i,int j)const {
  int min=MIN0(i,j);
  int max=MAX0(i,j);
  return(min+max*(max+1)/2);
}

double * basis_set::calc_olap() {
  int sh1, sh2;
  int n1, n2;
  int l1, l2;
  double * olap = new double[nbf * nbf];
  int shls[2];
  for (sh1 = 0; sh1 < nbas; sh1++) {
    n1 = 2 * bas[ANG_OF + BAS_SLOTS * sh1] + 1;
    for (sh2 = sh1; sh2 < nbas; sh2++) {
      n2 = 2 * bas[ANG_OF + BAS_SLOTS * sh2] + 1;
      double * buf = new double[n1 * n2];
      shls[0] = sh1;
      shls[1] = sh2;
      cint1e_ovlp_sph(buf, shls, atm, natm, bas, nbas, env);
      int i1, i2, j1, j2;
      for (i1 = first_bf_of_shell[sh1], j1 = 0; j1 < n1; i1++, j1++) {
	for (i2 = first_bf_of_shell[sh2], j2 = 0; j2 < n2; i2++, j2++) {
	  olap[i1 * nbf + i2] = buf[j1 * n2 + j2];
	  olap[i2 * nbf + i1] = buf[j1 * n2 + j2];
	}
      }
      delete[] buf;
    }
  }
  return (olap);
}
double * basis_set::calc_kin() {
  int sh1, sh2;
  int n1, n2;
  double * kin = new double[nbf * nbf];
  int shls[2];
  for (sh1 = 0; sh1 < nbas; sh1++) {
    int l1 = bas[ANG_OF + BAS_SLOTS * sh1];
    n1 = 2 * l1 + 1;
    for (sh2 = sh1; sh2 < nbas; sh2++) {
      int l2 = bas[ANG_OF + BAS_SLOTS * sh2];
      n2 = 2 * l2 + 1;
      double * buf = new double[n1 * n2];
      shls[0] = sh1;
      shls[1] = sh2;
      cint1e_kin_sph(buf, shls, atm, natm, bas, nbas, env);
      int i1, i2, j1, j2;
      for (i1 = first_bf_of_shell[sh1], j1 = 0; j1 < n1; i1++, j1++) {
	for (i2 = first_bf_of_shell[sh2], j2 = 0; j2 < n2; i2++, j2++) {
	  kin[i1 * nbf + i2] = buf[j1 * n2 + j2];
	  kin[i2 * nbf + i1] = buf[j1 * n2 + j2];
	}
      }
      delete[] buf;
    }
  }
  return (kin);
}
double * basis_set::calc_nuc() {
  int sh1, sh2;
  int n1, n2;
  double * nuc = new double[nbf * nbf];
  int shls[2];
  for (sh1 = 0; sh1 < nbas; sh1++) {
    int l1 = bas[ANG_OF + BAS_SLOTS * sh1];
    n1 = 2 * l1 + 1;
    for (sh2 = sh1; sh2 < nbas; sh2++) {
      int l2 = bas[ANG_OF + BAS_SLOTS * sh2];
      n2 = 2 * l2 + 1;
      double * buf = new double[n1 * n2];
      shls[0] = sh1;
      shls[1] = sh2;
      cint1e_nuc_sph(buf, shls, atm, natm, bas, nbas, env);
      int i1, i2, j1, j2;
      for (i1 = first_bf_of_shell[sh1], j1 = 0; j1 < n1; i1++, j1++) {
	for (i2 = first_bf_of_shell[sh2], j2 = 0; j2 < n2; i2++, j2++) {
	  nuc[i1 * nbf + i2] = buf[j1 * n2 + j2];
	  nuc[i2 * nbf + i1] = buf[j1 * n2 + j2];
	}
      }
      delete[] buf;
    }
  }
  return (nuc);
}
double * basis_set::calc_teint() {
  int sh1, sh2, sh3, sh4;
  int n1, n2, n3, n4;
  int nij = nbf * (nbf + 1) / 2;
  int nijkl = nij * (nij + 1) / 2;
  double * teint = new double[nijkl];
  int shls[4];
  for (sh1 = 0; sh1 < nbas; sh1++) {
    n1 = 2 * bas[ANG_OF + BAS_SLOTS * sh1] + 1;
    for (sh2 = sh1; sh2 < nbas; sh2++) {
      n2 = 2 * bas[ANG_OF + BAS_SLOTS * sh2] + 1;
      for (sh3 = 0; sh3 < nbas; sh3++) {
	n3 = 2 * bas[ANG_OF + BAS_SLOTS * sh3] + 1;
	for (sh4 = sh3; sh4 < nbas; sh4++) {
	  n4 = 2 * bas[ANG_OF + BAS_SLOTS * sh4] + 1;
		
	  // this prescreens condition ij<kl
	  int imax = first_bf_of_shell[sh1] + n1;
	  int jmax = first_bf_of_shell[sh2] + n2;
	  int ijmax = imax + jmax * (jmax + 1) / 2;
	  int kmin = first_bf_of_shell[sh3];
	  int lmin = first_bf_of_shell[sh4];
	  int klmin = kmin + lmin * (lmin + 1) / 2;
	  if (ijmax < klmin)
	    continue;
		
	  double * buf = new double[n1 * n2 * n3 * n4];
	  shls[0] = sh1;
	  shls[1] = sh2;
	  shls[2] = sh3;
	  shls[3] = sh4;
	  int ret = cint2e_sph(buf, shls, atm, natm, bas, nbas, env,
			       NULL);
	  if (ret == 0)
	    continue;
	  int i1, i2, i3, i4, j1, j2, j3, j4;
	  for (i1 = first_bf_of_shell[sh1], j1 = 0; j1 < n1;
	       i1++, j1++) {
	    for (i2 = first_bf_of_shell[sh2], j2 = 0; j2 < n2;
		 i2++, j2++) {
	      if (i1 > i2)
		continue;
	      int ij = packed_index(i1,i2);
	      for (i3 = first_bf_of_shell[sh3], j3 = 0; j3 < n3;
		   i3++, j3++) {
		for (i4 = first_bf_of_shell[sh4], j4 = 0;
		     j4 < n4; i4++, j4++) {
		  if (i3 > i4)
		    continue;
		  int kl = packed_index(i3,i4);
		  if (ij < kl)
		    continue;
		  int ijkl = packed_index(ij,kl);
		  teint[ijkl] = buf[j1 + j2 * n1
				    + j3 * n1 * n2 + j4 * n1 * n2 * n3];
		} //i4
	      } //i3
	    } //i2
	  } //i1
	  delete[] buf;
	} //sh4
      } //sh3
    } //sh2
  } //sh1
	
  return (teint);
}

double * basis_set::inv_sqrt(double * M) {
  double * V = new double[nbf * nbf];
  double * M2 = new double[nbf * nbf];

  double * eval = new double[nbf];
  memcpy(V, M, nbf * nbf * sizeof(double));
  int r = eigval(V, eval, nbf);


  for (int i = 0; i < nbf; i++) {
    for (int j = 0; j < nbf; j++) {
      M2[i * nbf + j] = 0.;
      for (int k = 0; k < nbf; k++) {
	M2[i * nbf + j] += V[k * nbf + i] * V[k * nbf + j]/ sqrt(eval[k]);
      }
    }
  }
  return (M2);
}

ostream& operator<<(ostream &out, const basis_set & b) {
  int i, j, n;
  out << "# Basis set information: " << endl;
  out << "# " << setw(30) << " number of atoms:" << setw(5) << b.natm << endl;

  out <<  "# Geometry in atomic units:" << endl;
  out <<  "# " << setw(5) <<"Z" << setw(15) << "x" << setw(15) << "y" << setw(15) << "z" << endl; 
    
  for (i = 0; i < b.natm; i++) {
    out <<  "# " << setw(5) << b.atm[CHARGE_OF + ATM_SLOTS * i]
	<< fixed << setw(15) << b.env[b.atm[PTR_COORD + ATM_SLOTS * i] + 0]
	<< fixed << setw(15) << b.env[b.atm[PTR_COORD + ATM_SLOTS * i] + 1]
	<< fixed << setw(15) << b.env[b.atm[PTR_COORD + ATM_SLOTS * i] + 2]
	<< endl; 
  }

  
  out << "# " << setw(30) << " number of basis functions:" << setw(5) << b.nbf << endl;  
  out << "# " << setw(30) << " number of shells:" << setw(5) << b.nbas << endl;

  out <<  "# Shells:" << endl;
  out <<  "# " << setw(10) <<"on atom" << setw(5) << "l" << setw(20) << "exp" << setw(20) << "coef" << endl; 

  for (n = 0; n < b.nbas; n++) {
    out <<  "# " << setw(10) << b.bas[ATOM_OF + BAS_SLOTS * n]
	<< setw(5) << b.bas[ANG_OF + BAS_SLOTS * n]
	<< setw(20) << ""
	<< setw(20) << ""
	<< endl; 
    for (j = 0; j < b.bas[NPRIM_OF + BAS_SLOTS * n]; j++) {
      out <<  "# " << setw(10) << ""
	  << setw(5) << ""
	  << setw(20) << b.env[b.bas[PTR_EXP + BAS_SLOTS * n] + j]
	  << setw(20) << b.env[b.bas[PTR_COEFF + BAS_SLOTS * n] + j]
	/CINTgto_norm(b.bas[ANG_OF + BAS_SLOTS * (n)],b.env[b.bas[PTR_EXP + BAS_SLOTS * n] + j])
	  << endl; 

    }
  }

  return (out);
}

