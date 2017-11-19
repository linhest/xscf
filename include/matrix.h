
int eigval(double * a, double * eval, int m);
void transform(double *a,double *b,int n);
double norm(double *C,double *S,int n);
void mmult(double *a,double *b,double * c,int n);
double trace(double * m,int n);
void print_matrix(double *m, int n, char * description);
int solve_AxB(double * a, double *b, double * d, int n);
