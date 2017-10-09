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

}
