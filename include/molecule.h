#include <vector>
#include <iostream>
using namespace std;

#define LEN_FILENAME 100

class molecule {
 public:
  molecule(char * filename);
  ~molecule();
  friend ostream& operator<< (ostream &out, const molecule &mol);
  char * get_basis_set();
  int get_natm();


  int natm; // Number of atoms
  double * geom; // geometry of atoms
  int * Z;// Z of atoms
  char ** atm_label;
  char basis_set[LEN_FILENAME];
  int nmo;// Number of orbitals
  int ms2;
  int nelec;
  int charge; // total charge

  double energy_conv;
  double density_conv;
  int n_diis; // dimension of diis space

  unsigned int * occ; // occupation of MOs

  double Enuc; // nuclear repulsion energy

  int print; // print level
  char filename_guess_save[LEN_FILENAME];
  char filename_guess_read[LEN_FILENAME];

 private:
  void read_basis();
  int parse_int(char * string);
  double parse_float(char * string);

  vector<char *>* parse_string_array(char * string, char delim);
  unsigned int * parse_int_array_fixed_len(unsigned int * s,  int nmo, char * buf, char delim);
  int * parse_int_array_fixed_len(int * s,  int nmo, char * buf, char delim);

  int nref; // number of reference occupations
  int exlevels; // number of excitaion levels
  int * endorbital;   // array of orbital indices: maximum orbital index per excitation level
  int * startorbital;  // array of orbital indices: minimal orbital index per excitation level 
  int * excludeorbital; // orbital indices to exclude in excitations
  int * fixorbital; // array of fixed occupations for each orbital
  int * rasorbital; // array: number of ras spaces x number of orbitals: 1 means orbital is insided ras space
  
  int ** rasocc; // ? 
  int nras; // number of ras spaces 
  int printoccs;
  unsigned int * targetocc; // nref x number of orbital indices : array containing reference occupations
  double spinl,spinm; // defines spin l and m 
  int target_sym; // symmetry irrep of solution to be found
  int numberofstates; // number of solutions to be calculated
  long int mem; // number of doubles to be allocated to store everything (if = -1 we look investigate how much we need automatically)
  unsigned int * refs;
};


static double calculate_Enuc(double * geom, int * Z, int natm);
