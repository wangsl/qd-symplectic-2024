
# $Id$

.DEFAULT: .F .F90 .c .C
.SUFFIXES: .F .F90 .c .C

SRC = $(PWD)
base_dir = $(SRC)/..
O = $(PWD)/../../build/HO2

F77 = gfortran
F90 = gfortran
CC = gcc
CXX = g++
NVCC = nvcc

CFLAGS = #-DI64

FFLAGS = -O3 -fPIC -fopenmp -fallow-argument-mismatch

CXXFLAGS = -O3 -fPIC -fopenmp -std=c++11 -I$(MATLAB_ROOT)/extern/include \
	-I$(base_dir)/common-cuda

Link = $(CXX)

LIBS = -L$(base_dir)/../build/common-cuda -lqmdyn -lgfortran

MEXA64Files = $(O)/DMBEIVMex.mexa64

OBJS = 	$(O)/dmbeiv.o $(O)/DMBEIVMex.o


.DEFAULT_GOAL := $(O)/DMBEIVMex.mexa64

all: $(MEXA64Files)

$(O)/%.o: %.c
	cd $(O) ; $(CC)  $(cFLAGS) -c $(SRC)/$<
$(O)/%.o: %.C
	cd $(O) ; $(CXX) $(CXXFLAGS) -c $(SRC)/$<
$(O)/%.o: %.F
	cd $(O) ; $(F77) $(FFLAGS) -c $(SRC)/$<
$(O)/%.o: %.F90
	cd $(O) ; $(F90) $(FFLAGS) -c $(SRC)/$<
$(O)/%.o: %.cu
	cd $(O) ; $(NVCC) $(NVCCFLAGS) -dc $(SRC)/$<

%io.C: %.h
	perl io.pl $<

$(O)/%.mexa64: $(OBJS)
	cd $(O); $(Link) -shared $(CXXFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o *~ *.mod $(EXE) depend $(MEXA64Files) $(QMLibs) $(OBJS)

cleancuda:
	rm -rf $(CUDAObjs) $(CUDALinkObj)

depend :
	$(CXX) $(CXXFLAGS) -MM *.[cC] | perl dep.pl > $@
	sh fortran-dep.sh *.F >> $@

ifneq ($(MAKECMDGOALS), clean)
include depend
endif

