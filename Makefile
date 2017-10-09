
DEBUG= -g -O0
#DEBUG=-g -O0 -check-pointers=rw -check-pointers-dangling=all -check-pointers-undimensioned #-lduma -DDUMA_PROTECT_BELOW
#DEBUG=  -p #-openmp-profile #-profile-functions #-profile-loops=all -profile-loops-report=2 

CC =  gcc
CPP = g++
FC = ifort2011

PSI3=/home/ludger/psi3
PSILDFLAGS=-L${PSI3}/lib/ -lPSI_chkpt -lPSI_iwl -lPSI_psio -lPSI_ciomr -lPSI_ipv1 -lPSI_qt 
#PSICFLAGS=-I${PSI3}/include/ -I${PSI3}/src/lib


#ARPACK= -L /home/ludger/ARPACK/ -larpack -lifcore -lifport -limf

#LDFLAGS = -O3 ${DEBUG} ${ARPACK} -openmp -openmp-report1 -lm -mkl=parallel

#CFLAGS = -O3 -openmp -openmp-report1 -I/home/linhest/include/ -mkl=parallel -O3 ${DEBUG} 

#TARGETS = csf_ci csf_ci_psi3 read_states read_csfs

#.PHONY: clean

molecule.o: molecule.cc
	${CPP} -g -c $< -o $@

basis_set.o: basis_set.cc
	${CPP} -g -c $< -o $@ -I/usr/local/include/


xscf.o: xscf.cc
	${CPP} -g -c $< -o $@ -I/usr/local/include/

xscf: xscf.o molecule.o basis_set.o 
	${CPP} -g  $^ -o $@  -L/usr/local/lib -lcint

all: ${TARGETS}

constants.o: constants.f90
	${FC} -c $< -o $@

anglib.o: anglib.f90
	${FC} -c $< -o $@

angularintegration.o: angularintegration.f90
	${FC} -c $< -o $@

readpsi3.o: readpsi3.cc types.h
	${CPP} ${CFLAGS} ${PSICFLAGS} -c  $< -o $@ 

%.o: %.c %.h types.h
	${CC} ${CFLAGS} -o $@ -c $<

csf_ci: csf_ci.o hamiltonian-spincsf.o occ.o det.o hamiltonian-det.o readmointegrals.o mointegrals.o  densitymatrix.o naturalorbitals.o
	${CPP}  -o $@ $^ ${LDFLAGS} 

csf_ci_psi3: csf_ci_psi3.o hamiltonian-spincsf.o occ.o det.o hamiltonian-det.o readpsi3.o mointegrals.o 
	${CPP}  -o $@ $^ ${LDFLAGS} ${PSILDFLAGS} 

print_naturalorbitals: print_naturalorbitals.o naturalorbitals.o densitymatrix.o 
	${CPP}  -o $@ $^ ${LDFLAGS} ${PSILDFLAGS} 

ionization_coef: ionization_coef.o 
	${CPP}  -o $@ $^ ${LDFLAGS} ${PSILDFLAGS} 

#auger_coef: auger_coef.o 
#	${CPP}  -o $@ $^ ${LDFLAGS} ${PSILDFLAGS} 

auger_coef_psi3: auger_coef_psi3.o readpsi3.o angularintegration.o constants.o anglib.o
	${CPP}  -o $@ $^ -lgfortran ${LDFLAGS} ${PSILDFLAGS} 

orbital_coef: orbital_coef.o readpsi3.o 
	${CPP}  -o $@ $^ ${LDFLAGS} ${PSILDFLAGS} 

overlap_states: overlap_states.o readpsi3.o 
	${CPP}  -o $@ $^ ${LDFLAGS} ${PSILDFLAGS} 

map_to_min_basis: map_to_min_basis.o readpsi3.o
	${CPP}  -o $@ $^ ${LDFLAGS} ${PSILDFLAGS} 

clean:	
	rm *.o 

.PHONY: debug


