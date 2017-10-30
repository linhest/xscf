int next_nhmp(int hold, int pold,int hold2, int pold2, int n,int m, unsigned int * occ,int orbend,int NORB,unsigned int * ref,unsigned int * list, int num);
unsigned int * excite(unsigned int * occ, int h, int p , int NORB );
unsigned int * add(unsigned int * occ, int p , int NORB );
unsigned int * rem(unsigned int * occ, int h , int NORB );
void print_occ(unsigned int *occ,int NORB);
void print_csf(spin_csf * csf );
void print_csf_2(spin_csf * csf ,double coef);

occ_set * occ_nhmp(molecule * mol,int n,int m,int orbstart, int orbend,unsigned int * ref);
spin_csf_set * occ2spincsf(unsigned int * occ, int NORB, double ls2, double ms,int * N);
det_set * occ2detset(unsigned int * occ,int NORB);

