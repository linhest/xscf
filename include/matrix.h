/*
 * matrix.h
 * written by Ludger Inhester
 * (c) 2017
 * 02.12.2017
 *
 */
int
eigval (double * a, double * eval, int m);
void
transform (double *a, double *b, int n);
double
norm (double *C, double *S, int n);
void
mmult (double *a, double *b, double * c, int n);
void
mmult_sym (double *a, double *b, double * c, int n);

double
trace (double * m, int n);
void
print_matrix (double *m, int n, char * description);
void
print_matrix (double *m, int n1, int n2, char * description);
int
solve_AxB (double * a, double *b, double * d, int n);
