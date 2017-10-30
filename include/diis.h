using namespace std;


class diis_set {
 public:
  diis_set(int nsave, int dim);
  ~diis_set();
 private:
  const int nsave;
  const int dim;
  //  vector <double*> Hsave; 
};
