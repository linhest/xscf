class basis_set{
public:
  basis_set(molecule & mol);
  friend ostream& operator<< (ostream &out, const basis_set &b);
  
private:
  int read_basis_set_atom_nshell(char * filename, char * atm_label);
  int read_basis_set_atom(char * filename, char * atm_label, double * env, int * bas, int atom_of, int *off, int *n);

  int natm;
  int nbas;
  int * atm;
  int * bas;
  double * env;
};

int symlabel2ang(char s);
