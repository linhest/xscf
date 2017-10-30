//determinant * neutral_determinant(molecule *mol);
det_set * build_det_set(int nelec,int casstart,int casend,molecule *mol);
det_set * add_electron(int i0, int nelec,int norb, det_set * dets);


void printstring(determinant * d);
int * addto_string(int * string, int o, int len );
int * copy_string(int * string, int len );
det_set * filter_det_set(det_set * dets);


int calc_det_sign(determinant * det,determinant * ref);
void calc_det_set_sign(det_set * dets);

det_set * remove_electron(determinant * ref);
determinant * groundstate_determinant(molecule *mol);
int * remove_from_string(int * string, int pos, int len );
