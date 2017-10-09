

typedef struct{
  int nmo;// Number of orbitals
  int ms2;
  int nelec;
  int charge; // total charge
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
}molecule;



typedef struct{
  int nhole;// Number of holes
  int npart;// Number of particle
  int * hole;
  int * part;
}configuration_state_function;

typedef struct{
  int n;
  configuration_state_function ** csf;
}csf_set;

typedef struct{
  int n;
  unsigned int * occ;
}occ_set;

typedef struct{
  //  int NORB;
  //int NELEC;
  int nalpha;
  int nbeta;
  int * alphastring;
  int * betastring;
}determinant;

typedef struct{
  int n;
  int NORB;
  determinant ** det;
  //  determinant * ref;
}det_set;

typedef struct{
  //  unsigned int degeneracy;
  //  double l2;
  //  double m;
  int noccs;
  unsigned int * occ;
  //  int n;
  det_set * dets;
  double * c;
}spin_csf;

typedef struct{
  int n;
  spin_csf ** csf;
}spin_csf_set;






