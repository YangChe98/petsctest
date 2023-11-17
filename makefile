FILES = Makefile main.c 
#PETSC_DIR = /usr/lib/petscdir/petsc64-3.15 
default: all

all: main

clean:
	-/bin/rm -f *.o *.tmp plot.eps PI[0-9]* \
                *.G.* *.A.* *.Aalpha.* *.Abeta.* \
                *.Gx.* *.Gy.* *.Gz.* *.x.* *.y.* *.z.* *.b.* *.x0.*

veryclean: clean
	-/bin/rm -f poisson simplest maxwell heat mesh*.* eigen elastic \
                maxwell-complex maxwell-eigen maxwell-eigen1 non-smooth \
                navier-stokes main *.vtk

#include /share/soft/phg/phg-0.9.4-mvapich2-20190318/share/phg/Makefile.inc
#include $(PHG_MAKEFILE_INC)
#include /soft/apps/phg/gcc-10.2.0/mvapich2-2.3.5/phg-0.9.6/share/phg/Makefile.inc
#include /usr/lib/petscdir/petsc64-3.15/x86_64-linux-gnu-real/include/petscksp.h 
#include /usr/lib/petscdir/petsc64-3.15/x86_64-linux-gnu-real/share/petsc/Makefile.user

#main.o: main.c

main: main.c
	mpicc -o main main.c -I/usr/lib/petscdir/petsc64-3.15/x86_64-linux-gnu-real/include/  -lpetsc
#main.o: main.c

#main: main.o
#	${LINKER} ${LDFLAGS} -o $@ $^ ${LIBS}
.PHONY: default all clean veryclean lib
                                           