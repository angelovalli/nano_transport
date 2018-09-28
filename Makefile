##$ COMPILER: suppprted compilers are ifort, gnu >v4.7 or use mpif90
FC=gfortran

##$ PLATFORM: supported platform are intel, gnu 
PLAT=gnu

##$ SET THE LOCATION OF YOU PROGRAM DRIVER (default is ./drivers)
DIR =.

##$ SET THE TARGET DIRECTORY WHERE TO PUT THE EXECUTABLE (default if $HOME/.bin in the PATH)
DIREXE=$(HOME)/.bin

##$ CHOOSE THE DRIVER CODE:
EXE=main

##$ SET INCLUDE AND LINK OPTIONS USING pkg-config
# INCARGS=$(shell pkg-config --cflags blas lapack)
# LIBARGS=$(shell pkg-config --libs   blas lapack)


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



###############################################################################
#                        HERE STARTS THE REAL WORK
###############################################################################
OBJS=INPUT_VARS.o COMMON.o BUILD_H.o GREENS_FUNCTIONS.o TRANSPORT.o
BLASOBJS=dcabs1.o ztrmm.o ztrmv.o izamax.o lsame.o  xerbla.o zgemm.o  zgemv.o  zscal.o  zswap.o  ztrsm.o
LAPACKOBJS=dlamch.o  ieeeck.o  ilaenv.o  iparmq.o  zgetrf.o  zgetrf2.o zgetri.o  zlaswp.o  ztrti2.o  ztrtri.o

all: lib compile

debug: debug lib_debug compile
debug: FFLAG=$(DFLAG)


compile: $(OBJS)
	@echo " ..................... compile ........................... "
	$(FC) $(FFLAG) $(INCARGS) $(OBJS) $(BLASOBJS) $(LAPACKOBJS) $(DIR)/$(EXE).f90 -o $(DIREXE)/transport_$(EXE)
	@echo " ...................... done .............................. "
	@echo ""
	@echo ""
	@echo "created" $(DIREXE)/transport_$(EXE)

lib: 
	@make -C blas
	@make -C lapack



lib_debug:
	@make -C blas debug 
	@make -C lapack debug


.f90.o:	
	$(FC) $(FFLAG) $(INCARGS) -c $<


clean: 
	@echo "Cleaning"
	@rm -f *.mod *.o *~

allclean: clean
	@make -C blas clean
	@make -C lapack clean


