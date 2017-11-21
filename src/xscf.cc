#include <iostream>
#include <deque>

#include "xscf.h"
#include "matrix.h"
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
  cout << mol;
  basis_set basis(mol);
  cout << basis << endl;

  if(mol.print>2){
      print_matrix(basis.olap,basis.nbf,(char *)"Overlap Matrix");
      print_matrix(basis.olap_inv_sqrt,basis.nbf,(char *)"Overlap Matrix inv sqrt");
      print_matrix(basis.kin,basis.nbf,(char *)"Hkin Matrix");
      print_matrix(basis.nuc,basis.nbf,(char *)"nuc Matrix");
  }

  scf S(mol,basis);
  if(mol.filename_guess_save[0]!='\0')
    S.save_guess(mol.filename_guess_save);


}
