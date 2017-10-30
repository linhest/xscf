

class scf {
 public:
  scf(const molecule &arg_mol, const basis_set &arg_b);
  ~scf();

 private:
  int print_iteration();
  double total_energy(double * dmat,double * h0, double * f) const;
  int construct_dmat(double  * dmat, double * c, const unsigned int * occ );
  void construct_fmat(double * f, double * h0, double * dmat ,double * teint);
  const molecule & mol;
  const basis_set & b;
  diis_set diis;
  double energy;
  double old_energy;


  double * c;
  double * dmat;

  double * oe;
  double * h0;
  double * f;
  int it;
};


