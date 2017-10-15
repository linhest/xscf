
DEBUG= -g -O0
#DEBUG=-g -O0 -check-pointers=rw -check-pointers-dangling=all -check-pointers-undimensioned #-lduma -DDUMA_PROTECT_BELOW
#DEBUG=  -p #-openmp-profile #-profile-functions #-profile-loops=all -profile-loops-report=2 

CC =  gcc
CPP = g++
FC = ifort2011

#PSI3=/home/ludger/psi3
#PSILDFLAGS=-L${PSI3}/lib/ -lPSI_chkpt -lPSI_iwl -lPSI_psio -lPSI_ciomr -lPSI_ipv1 -lPSI_qt 
#PSICFLAGS=-I${PSI3}/include/ -I${PSI3}/src/lib


#ARPACK= -L /home/ludger/ARPACK/ -larpack -lifcore -lifport -limf

#LDFLAGS = -O3 ${DEBUG} ${ARPACK} -openmp -openmp-report1 -lm -mkl=parallel

#CFLAGS = -O3 -openmp -openmp-report1 -I/home/linhest/include/ -mkl=parallel -O3 ${DEBUG} 

TARGETS = xscf

all: ${TARGETS}


.PHONY: clean

matrix.o: matrix.cc 
	${CPP} -g -c $< -o $@

molecule.o: molecule.cc
	${CPP} -g -c $< -o $@

scf.o: scf.cc
	${CPP} -g -c $< -o $@


basis_set.o: basis_set.cc
	${CPP} -g -c $< -o $@ -I/usr/local/include/ -I/usr/include/


xscf.o: xscf.cc
	${CPP} -g -c $< -o $@ -I/usr/local/include/

xscf: xscf.o molecule.o basis_set.o matrix.o scf.o
	${CPP} -g  $^ -o $@  -L/usr/local/lib -lcint -llapack -lblas



clean:	
	rm *.o 

.PHONY: debug


