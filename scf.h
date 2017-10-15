

class scf {
 public:
	scf(const molecule &arg_mol, const basis_set &arg_b);
  ~scf();

 private:
  int print_iteration();
  int construct_dmat(double  * dmat, double * c, const unsigned int * occ );
  const molecule & mol;
  const basis_set & b;
  double energy;
  double * c;
  double * dmat;

  double * oe;
  double * h0;
  double * f;
};


