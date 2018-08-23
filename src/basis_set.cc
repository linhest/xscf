/*
 * basis_set.cc
 * written by Ludger Inhester
 * (c) 2017
 * 02.12.2017
 *
 */
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include "../include/molecule.h"

#include "../include/basis_set.h"
#include "../include/matrix.h"

#include <cstdlib>
using namespace std;

extern "C" {
#include <cint.h>
  void
  CINTinit_2e_optimizer (CINTOpt **opt, const FINT *atm, const FINT natm,
                         const FINT *bas, const FINT nbas, const double *env);

  int
  cint1e_ovlp_sph (double *buf, int *shls, int *atm, int natm, int *bas,
                   int nbas, double *env);
  int
  cint1e_kin_sph (double *buf, int *shls, int *atm, int natm, int *bas,
                  int nbas, double *env);
  int
  cint1e_nuc_sph (double *buf, int *shls, int *atm, int natm, int *bas,
                  int nbas, double *env);

  int
  dgemm_ (char *transa, char *transb, int *m, int * n, int *k, double *alpha,
          double *a, int *lda, double *b, int *ldb, double *beta, double *c,
          int *ldc);
  int
  dgemv_ (char *trans, int *m, int *n, double * alpha, double *a, int *lda,
          double *x, int *incx, double *beta, double *y, int *incy);
}
#define ENV_SIZE 10000
#define MAXLINELENGTH 2000

basis_set::~basis_set () {
  delete[] atm;
  delete[] bas;
  delete[] env;
  delete[] first_bf_of_shell;
  delete[] shells_on_atom;
  delete[] olap;
  delete[] kin;
  delete[] nuc;
  delete[] teint;
  delete[] olap_inv_sqrt;
}

// initializer of basis_set object
// reads in all basis_set parameters (exponents and contraction coefficients)
basis_set::basis_set (const molecule & mol) : mol (mol) {
  const char * filename = &mol.basis_set[0] ;
  int off = PTR_ENV_START; // = 20
  int n = 0;
  for (int i = 0; i < mol.natm; i++) {
    read_basis_set_atom (filename, mol.atm_label[i], NULL, NULL, i, &off, &n);
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
    read_basis_set_atom (filename, mol.atm_label[i], env, bas, i, &off, &n);
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
    first_bf_of_shell[n] = nbf;
    nbf += 2 * l + 1;
  }

  /*  if (nbf != mol.nmo) {
    char text[100];
    sprintf (text, "nmo should be %d", nbf);
    throw ::std::invalid_argument (text);
  }
  */

  // normalize basis functions
  if (bas != NULL) {
    int sh1;
    int n1;
    int shls[2];
    for (sh1 = 0; sh1 < nbas; sh1++) {
      n1 = 2 * bas[ANG_OF + BAS_SLOTS * sh1] + 1;
      double * buf = new double[n1 * n1];
      shls[0] = sh1;
      shls[1] = sh1;
      cint1e_ovlp_sph (buf, shls, atm, natm, bas, nbas, env);
      int i1, j1;
      if (mol.print > 1) {
        for (i1 = first_bf_of_shell[sh1], j1 = 0; j1 < n1; i1++, j1++) {
          printf ("shell %d bf %d norm %e\n", sh1, i1, buf[j1 * n1 + j1]);
        }
      }
      for (int e = 0; e < bas[NPRIM_OF + BAS_SLOTS * sh1]; e++) {
        env[bas[PTR_COEFF + BAS_SLOTS * sh1] + e] /= sqrt (buf[0]);
      }
      delete[] buf;
    }
  }
  // calculate quantities for this basis set
  olap = calc_olap ();
  kin = calc_kin ();
  nuc = calc_nuc ();
  teint = calc_teint ();
  olap_inv_sqrt = inv_sqrt (olap);
}

/* read_basis_set_atom
 reads in basis set data
 from file filename
 for atom atm_label with index atom_of
 env,bas,off,n are the libcint arrays
 */
int
basis_set::read_basis_set_atom (const char * filename, char * atm_label, double * env,
                                int * bas, int atom_of, int *off, int *n) {
  ifstream fin;
  bool found = false;
  int nshell_atm = 0;
  fin.open (filename);
  if (!fin.good ())
    throw std::invalid_argument (" basis filename not found ");
  while (!fin.eof ()) {
    char alabel[4];
    char slabel[2];
    int nexp;
    double exp, norm, coef, pcoef;
    int ai;
    int num;
    char buf[MAXLINELENGTH];
    fin.getline (buf, MAXLINELENGTH);
    if (strncmp ("****", buf, 4) == 0) {
      // new entry starts here
      fin.getline (buf, MAXLINELENGTH);
      if (fin.eof ())
        break;
      num = sscanf (buf, "%s %d", alabel, &ai);
      if (num != 2)
        throw std::invalid_argument (strcat (buf, "alabel int"));
    } else {
      num = sscanf (buf, "%s %d %lf", slabel, &nexp, &norm);
      if (num != 3)
        throw std::invalid_argument (strcat (buf, "sym nexp norm"));
      if (strcasecmp (alabel, atm_label) == 0) {
        found = true;
        if (strlen (slabel) == 1) {
          nshell_atm++;
          if (bas != NULL && env != NULL) {
            bas[ATOM_OF + BAS_SLOTS * *n] = atom_of;
            bas[ANG_OF + BAS_SLOTS * *n] = symlabel2ang (slabel[0]);
            bas[NPRIM_OF + BAS_SLOTS * *n] = nexp;
            bas[NCTR_OF + BAS_SLOTS * *n] = 1;
            bas[PTR_EXP + BAS_SLOTS * *n] = *off;
            bas[PTR_COEFF + BAS_SLOTS * *n] = *off + nexp;
          }
          for (int i = 0; i < nexp; i++) {
            fin.getline (buf, MAXLINELENGTH);
            num = sscanf (buf, "%lf %lf ", &exp, &coef);
            if (num != 2)
              throw std::invalid_argument (strcat (buf, "exp coef"));
            if (bas != NULL && env != NULL) {
              env[*off] = exp;
              env[*off + nexp] = coef
                  * CINTgto_norm (bas[ANG_OF + BAS_SLOTS * (*n)], exp);
            }
            (*off)++;
          }
          (*off) += nexp;
          (*n)++;
        } else if (strlen (slabel) == 2 && strcmp (slabel, "SP") == 0) {
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
            fin.getline (buf, MAXLINELENGTH);
            num = sscanf (buf, "%lf %lf %lf", &exp, &coef, &pcoef);
            if (num != 3)
              throw std::invalid_argument (strcat (buf, "exp scoef pcoef"));
            if (bas != NULL && env != NULL) {
              env[*off] = exp;
              env[*off + nexp] = coef * CINTgto_norm (0, exp);
              env[*off + 2 * nexp] = pcoef * CINTgto_norm (1, exp);
            }
            (*off)++;
          }
          (*off) += 2 * nexp;
        }
      } else {
        if (found)
          break;
        for (int i = 0; i < nexp; i++) {
          fin.getline (buf, MAXLINELENGTH);
        }
      }
    }
  }
  fin.close ();
  return (*n);
}

/*
 * converts "SPDFGHI" into 01234567
 */
int
symlabel2ang (char s) {
  char slabel[] = "SPDFGHI";
  for (int i = 0; i < 7; i++) {
    if (slabel[i] == s)
      return (i);
  }
  throw std::invalid_argument (" unknown angular momentum ");
  return(-1);
}

#define MIN0(a,b) (((a)<(b)) ? (a) : (b))
#define MAX0(a,b) (((a)>(b)) ? (a) : (b))
// index for packed storage of symmetric matrices
int
basis_set::packed_index (int i, int j) const {
  int min = MIN0(i, j);
  int max = MAX0(i, j);
  return (min + max * (max + 1) / 2);
}

/*
 *  calculates overlap matrix S=< phi_i | phi_j >
 */
double *
basis_set::calc_olap () {
  int sh1, sh2;
  int n1, n2;
  double * olap = new double[nbf * nbf];
  memset(olap, 0, nbf*nbf*sizeof(double));
  int shls[2];
  for (sh1 = 0; sh1 < nbas; sh1++) {
    n1 = 2 * bas[ANG_OF + BAS_SLOTS * sh1] + 1;
    for (sh2 = sh1; sh2 < nbas; sh2++) {
      n2 = 2 * bas[ANG_OF + BAS_SLOTS * sh2] + 1;
      double * buf = new double[n1 * n2];
      shls[0] = sh1;
      shls[1] = sh2;
      int ret = cint1e_ovlp_sph (buf, shls, atm, natm, bas, nbas, env);
      if (ret == 0)
        continue;
      int i1, i2, j1, j2;
      for (i2 = first_bf_of_shell[sh2], j2 = 0; j2 < n2; i2++, j2++) {
        for (i1 = first_bf_of_shell[sh1], j1 = 0; j1 < n1; i1++, j1++) {
          olap[i1 * nbf + i2] = buf[j1 + n1 * j2];
          olap[i2 * nbf + i1] = buf[j1 + n1 * j2];
        }
      }
      delete[] buf;
    }
  }
  return (olap);
}
/*
 *  calculates < phi_i | hkin | phi_j >
 */
double *
basis_set::calc_kin () {
  int sh1, sh2;
  int n1, n2;
  double * kin = new double[nbf * nbf];
  memset(kin, 0, nbf*nbf*sizeof(double));
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
      int ret = cint1e_kin_sph (buf, shls, atm, natm, bas, nbas, env);
      if (ret == 0)
        continue;
      int i1, i2, j1, j2;
      for (i2 = first_bf_of_shell[sh2], j2 = 0; j2 < n2; i2++, j2++) {
        for (i1 = first_bf_of_shell[sh1], j1 = 0; j1 < n1; i1++, j1++) {
          kin[i1 * nbf + i2] = buf[j1 + n1 * j2];
          kin[i2 * nbf + i1] = buf[j1 + n1 * j2];
        }
      }
      delete[] buf;
    }
  }
  return (kin);
}
/*
 *  calculates < phi_i | V_nuc | phi_j >
 */
double *
basis_set::calc_nuc () {
  int sh1, sh2;
  int n1, n2;
  double * nuc = new double[nbf * nbf];
  memset(nuc, 0, nbf*nbf*sizeof(double));
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
      int ret = cint1e_nuc_sph (buf, shls, atm, natm, bas, nbas, env);
      if (ret == 0)
        continue;
      int i1, i2, j1, j2;
      for (i2 = first_bf_of_shell[sh2], j2 = 0; j2 < n2; i2++, j2++) {
        for (i1 = first_bf_of_shell[sh1], j1 = 0; j1 < n1; i1++, j1++) {
          nuc[i1 * nbf + i2] = buf[j1 + n1 * j2];
          nuc[i2 * nbf + i1] = buf[j1 + n1 * j2];
        }
      }
      delete[] buf;
    }
  }
  return (nuc);
}
/*
 *  calculates all two electron integrals
 */
double *
basis_set::calc_teint () {
  int sh1, sh2, sh3, sh4;
  int n1, n2, n3, n4;
  int nij = nbf * (nbf + 1) / 2;
  int nijkl = nij * (nij + 1) / 2;
  double * teint = new double[nijkl];
  memset(teint, 0, nijkl*sizeof(double));
  int shls[4];
  CINTOpt *opt = NULL;
  cint2e_sph_optimizer (&opt, atm, natm, bas, nbas, env);
  for (sh1 = 0; sh1 < nbas; sh1++) {
    n1 = 2 * bas[ANG_OF + BAS_SLOTS * sh1] + 1;
    for (sh2 = 0; sh2 <= sh1; sh2++) {
      n2 = 2 * bas[ANG_OF + BAS_SLOTS * sh2] + 1;
      for (sh3 = 0; sh3 < nbas; sh3++) {
        n3 = 2 * bas[ANG_OF + BAS_SLOTS * sh3] + 1;
        for (sh4 = 0; sh4 <= sh3; sh4++) {
          n4 = 2 * bas[ANG_OF + BAS_SLOTS * sh4] + 1;
          // this prescreens condition ij<kl
          int imax = first_bf_of_shell[sh1] + n1;
          int jmax = first_bf_of_shell[sh2] + n2;
          int ijmax = jmax + imax * (imax + 1) / 2;
          int kmin = first_bf_of_shell[sh3];
          int lmin = first_bf_of_shell[sh4];
          int klmin = lmin + kmin * (kmin + 1) / 2;
          if (ijmax < klmin)
            continue;
          double * buf = new double[n1 * n2 * n3 * n4];
          shls[0] = sh1;
          shls[1] = sh2;
          shls[2] = sh3;
          shls[3] = sh4;
          int ret = cint2e_sph (buf, shls, atm, natm, bas, nbas, env, opt);
          if (ret == 0)
            continue;
          int i1, i2, i3, i4, j1, j2, j3, j4;
          for (i4 = first_bf_of_shell[sh4], j4 = 0; j4 < n4; i4++, j4++) {
            for (i3 = first_bf_of_shell[sh3], j3 = 0; j3 < n3; i3++, j3++) {
              if (i3 < i4)
                continue;
              int kl = packed_index (i3, i4);
              for (i2 = first_bf_of_shell[sh2], j2 = 0; j2 < n2; i2++, j2++) {
                for (i1 = first_bf_of_shell[sh1], j1 = 0; j1 < n1; i1++, j1++) {
                  if (i1 < i2)
                    continue;
                  int ij = packed_index (i1, i2);
                  if (ij < kl)
                    continue;
                  int ijkl = packed_index (ij, kl);
                  //                  printf("(%d %d|%d %d)] %d\n",i1,i2,i3,i4,ijkl);

                  teint[ijkl] = buf[j1 + j2 * n1 + j3 * n1 * n2
                                    + j4 * n1 * n2 * n3];
                } //i1
              } //i2
            } //i3
          } //i4
          delete[] buf;
        } //sh4
      } //sh3
    } //sh2
  } //sh1
  CINTdel_optimizer (&opt);
  return (teint);
}

/*
 *  calculates M=S^{-1/2}
 */
double *
basis_set::inv_sqrt (double * M) {
  double * V = new double[nbf * nbf];
  double * M2 = new double[nbf * nbf];

  double * eval = new double[nbf];
  memcpy (V, M, nbf * nbf * sizeof(double));
  eigval (V, eval, nbf);
  printf ("# smallest S eigenvalue %e\n", eval[0]);
  for (int i = 0; i < nbf; i++) {
    for (int j = 0; j < nbf; j++) {
      M2[i * nbf + j] = 0.;
      for (int k = 0; k < nbf; k++) {
        M2[i * nbf + j] += V[k * nbf + i] * V[k * nbf + j] / sqrt (eval[k]);
      }
    }
  }
  delete[] eval;
  delete[] V;
  return (M2);
}

/*
 * print out basis set information
 */
ostream&
operator<< (ostream &out, const basis_set & b) {
  int i, j, n;
  out << "# Basis set information: " << endl;
  out << "# " << setw (30) << " number of atoms:" << setw (5) << b.natm << endl;
  out << "# Geometry in atomic units:" << endl;
  out << "# " << setw (5) << "Z" << setw (15) << "x" << setw (15) << "y"
      << setw (15) << "z" << endl;
  for (i = 0; i < b.natm; i++) {
    out << "# " << setw (5) << b.atm[CHARGE_OF + ATM_SLOTS * i] << fixed
        << setw (15) << b.env[b.atm[PTR_COORD + ATM_SLOTS * i] + 0] << fixed
        << setw (15) << b.env[b.atm[PTR_COORD + ATM_SLOTS * i] + 1] << fixed
        << setw (15) << b.env[b.atm[PTR_COORD + ATM_SLOTS * i] + 2] << endl;
  }
  out << "# " << setw (30) << " number of basis functions:" << setw (5) << b.nbf
      << endl;
  out << "# " << setw (30) << " number of shells:" << setw (5) << b.nbas
      << endl;
  if (b.mol.print > 1) {
    out << "# Shells:" << endl;
    out << "# " << setw (10) << "on atom" << setw (5) << "l" << setw (20)
                                    << "exp" << setw (20) << "coef" << endl;
    for (n = 0; n < b.nbas; n++) {
      out << "# " << setw (10) << b.bas[ATOM_OF + BAS_SLOTS * n] << setw (5)
                                      << b.bas[ANG_OF + BAS_SLOTS * n] << setw (20) << "" << setw (20) << ""
                                      << endl;
      for (j = 0; j < b.bas[NPRIM_OF + BAS_SLOTS * n]; j++) {
        out << "# " << setw (10) << "" << setw (5) << "" << setw (20)
                                        << b.env[b.bas[PTR_EXP + BAS_SLOTS * n] + j] << setw (20)
                                        << b.env[b.bas[PTR_COEFF + BAS_SLOTS * n] + j]
                                                 / CINTgto_norm (b.bas[ANG_OF + BAS_SLOTS * (n)],
                                                                 b.env[b.bas[PTR_EXP + BAS_SLOTS * n] + j])
                                                                 << endl;
      }
    }
  }
  return (out);
}

