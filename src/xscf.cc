#include <iostream>
#include "molecule.h"
#include "basis_set.h"
#include "diis.h"
#include "scf.h"

using namespace std;
int main(int argc, char ** argv){
	char * input_filename;
	char  default_input_filename[]= "input.txt";
	if(argc==1){
	//	 throw std::invalid_argument( "please provide input file!" );
		cout << default_input_filename << endl;
		input_filename=default_input_filename;
	}else{
		input_filename=argv[1];

	}
  molecule mol(input_filename);
  cout << "test" << endl;
  cout << mol;
  basis_set basis(mol);
  cout << basis << endl;

  /*  cout << "#Overlap Matrix:" << endl;
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

  cout << "end " << endl;
  */

  scf S(mol,basis);

}
