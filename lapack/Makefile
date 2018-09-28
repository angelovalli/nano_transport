##$ COMPILER: suppprted compilers are ifort, gnu >v4.7 or use mpif90
FC=gfortran

##$ PLATFORM: supported platform are intel, gnu 
PLAT=gnu

ifeq ($(PLAT),intel)
FFLAG=-O2 -ftz
OFLAG=-O3 -ftz
DFLAG=-p -O0 -g -fpe0 -warn -warn errors -debuEg extended -traceback -check all,noarg_temp_created
endif

ifeq ($(PLAT),gnu)
FFLAG = -O2 -ffree-line-length-none
DFLAG = -O0 -p -g -fimplicit-none -Wsurprising  -Waliasing -fwhole-file -fcheck=all -pedantic -fbacktrace -ffree-line-length-none
OFLAG = -O3 -ffast-math -march=native -funroll-all-loops -fno-protect-parens -flto -ffree-line-length-none
endif

##$ Extends the implicit support of the Makefile to .f90 files
.SUFFIXES: .f90 .f


LAPACKLIB=dlamch.o  ieeeck.o  ilaenv.o  iparmq.o  zgetrf.o  zgetrf2.o zgetri.o  zlaswp.o  ztrti2.o  ztrtri.o

all: lib

debug: debug lib
debug: FFLAG=$(DFLAG)

lib: lapack

lapack: $(LAPACKLIB)
	@echo "LAPACK done"
	@rsync -a *.o ../
	@echo ""

.f.o:	
	$(FC) $(FFLAG) $(INCARGS) -c $<


clean: 
	@echo "Cleaning LAPACK:"
	@rm -fv *.o


