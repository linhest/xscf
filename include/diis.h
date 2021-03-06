/*
 * diis.h
 * written by Ludger Inhester
 * (c) 2017
 * 02.12.2017
 *
 */
using namespace std;
#include <deque>

class diis_set {
public:
  diis_set (int nsave, int dim);
  double *
  diis_interpolate (double * dmat, double * fmat);
  double
  diis_error (double * dmat, double * fmat);

  ~diis_set ();
private:
  //  double * olap;
  const int nsave;
  const int dim;
  deque<double*> fSave;
  deque<double*> errSave;

};

