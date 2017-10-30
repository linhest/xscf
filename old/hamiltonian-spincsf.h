double hamiltonian_spincsf_entry(spin_csf * left,
                                 spin_csf * right);
void build_hamiltonian_spincsf(double *h,spin_csf_set * scsf);

double * diagonal_hamiltonian_spincsf(spin_csf_set *csfs);
long int build_sparse_hamiltonian_spincsf(spin_csf_set * scsf,double * val, int * rowind, int * colind);
