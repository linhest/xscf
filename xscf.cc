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

  cout << "#Overlap Matrix:" << endl;
  for(int i=0;i<basis.nbf;i++){
	  for(int j=0;j<basis.nbf;j++){
		  cout << basis.olap[i*basis.nbf+j]<< " ";
	  }
	  cout << endl;
  }

  cout << "#Overlap Matrix inv sqrt:" << endl;
  for(int i=0;i<basis.nbf;i++){
	  for(int j=0;j<basis.nbf;j++){
		  cout << basis.olap_inv_sqrt[i*basis.nbf+j]<< " ";
	  }
	  cout << endl;
  }

  cout << "#Hkin Matrix:" << endl;
  for(int i=0;i<basis.nbf;i++){
	  for(int j=0;j<basis.nbf;j++){
		  cout << basis.kin[i*basis.nbf+j]<< " ";
	  }
	  cout << endl;
  }


  cout << "#nuc Matrix:" << endl;
  for(int i=0;i<basis.nbf;i++){
	  for(int j=0;j<basis.nbf;j++){
		  cout << basis.nuc[i*basis.nbf+j]<< " ";
	  }
	  cout << endl;
  }


  double * h0 =new double[basis.nbf*basis.nbf];
  for(int i=0;i<basis.nbf*basis.nbf;i++){
	  h0[i]=basis.kin[i]+basis.nuc[i];
  }

  delete[] h0;

}
