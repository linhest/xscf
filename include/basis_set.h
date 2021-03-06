/*
 * basis_set.h
 * written by Ludger Inhester
 * (c) 2017
 * 02.12.2017
 *
 */
class basis_set {
public:
  basis_set (const molecule & mol);
  ~basis_set ();
  friend ostream&
  operator<< (ostream &out, const basis_set &b);
  int
  packed_index (int i, int j) const;

  const molecule & mol;

  double * teint;
  double * kin;
  double * olap;
  double * nuc;
  double * olap_inv_sqrt;

  int natm;
  int nbas;
  int * atm;
  int * bas;
  double * env;
  int * shells_on_atom;
  int * first_bf_of_shell;
  int nbf;
  int lmax;
private:
  int
  read_basis_set_atom (const char * filename, char * atm_label, double * env,
                       int * bas, int atom_of, int *off, int *n);

  double *
  calc_olap ();
  double *
  calc_kin ();
  double *
  calc_nuc ();
  double *
  calc_teint ();
  double *
  inv_sqrt (double * M);

};

int
symlabel2ang (char s);
