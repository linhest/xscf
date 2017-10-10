#include <iostream>
#include "molecule.h"
#include "basis_set.h"

using namespace std;
int main(int argc, char ** argv){
  molecule mol(argv[1]);
  cout << "test" << endl;
  cout << mol;
  basis_set basis(mol);
  cout << basis << endl;

  double * olap=basis.olap();
  cout << "#Overlap Matrix:" << endl;
  for(int i=0;i<basis.nbf;i++){
	  for(int j=0;j<basis.nbf;j++){
		  cout << olap[i*basis.nbf+j]<< " ";
	  }
	  cout << endl;
  }
  delete[] olap;


  double * kin=basis.kin();
  cout << "#Hkin Matrix:" << endl;
  for(int i=0;i<basis.nbf;i++){
	  for(int j=0;j<basis.nbf;j++){
		  cout << kin[i*basis.nbf+j]<< " ";
	  }
	  cout << endl;
  }
  delete[] kin;


  double * nuc=basis.nuc();
  cout << "#nuc Matrix:" << endl;
  for(int i=0;i<basis.nbf;i++){
	  for(int j=0;j<basis.nbf;j++){
		  cout << nuc[i*basis.nbf+j]<< " ";
	  }
	  cout << endl;
  }
  delete[] nuc;

}
