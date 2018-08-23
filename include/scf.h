/*
 * scf.h
 * written by Ludger Inhester
 * (c) 2017
 * 02.12.2017
 *
 */

class scf {
public:
  scf (const molecule &arg_mol, const basis_set &arg_b);
  ~scf ();
  int
  save_guess (const char * filename);
private:
  void
  print_iteration ();
  void
  print_result ();

  int
  read_guess (const char * filename);

  double
  total_energy (double * dmat, double * h0, double * f) const;
  double
  diis_error (double * dmat, double *f, double * err);
  void
  construct_dmat (double * dmat, double * c, const unsigned int * occ);
  void
  construct_fmat (double * f, double * h0, double * dmat, double * teint);
  const molecule & mol;
  const basis_set & b;
  diis_set diis;
  double energy;
  double old_energy;
  double max_dmat_diff;

  double * c;
  double * dmat;
  double * dmat_old;
  unsigned int * occ;
  double * oe;
  double * h0;
  double * f;
  int it;
};

