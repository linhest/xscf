class basis_set{
public:
  basis_set(molecule & mol);
  ~basis_set();
  friend ostream& operator<< (ostream &out, const basis_set &b);
  double * olap();
  double * kin();
  double * nuc();
  double * teint();

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
  int read_basis_set_atom_nshell(char * filename, char * atm_label);
  int read_basis_set_atom(char * filename, char * atm_label, double * env, int * bas, int atom_of, int *off, int *n);


};

int symlabel2ang(char s);
