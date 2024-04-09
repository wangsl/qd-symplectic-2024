
# $Id$

.DEFAULT: .F .F90 .c .C
.SUFFIXES: .F .F90 .c .C

SRC = $(PWD)
O = $(SRC)/../../build/H3

F77 = gfortran
F90 = gfotran
CC = gcc
CXX = g++

CFLAGS = -O3 -fopenmp -fPIC 

CXXFLAGS = -O3 -fopenmp -fPIC -std=c++0x -I$(MATLAB_ROOT)/extern/include -I$(SRC)/../common-cuda

FFLAGS = -O3 -fPIC -fopenmp -ffixed-line-length-132

Link = $(CXX)

LIBS =  -L$(SRC)/../../build/common-cuda -lqmdyn -lgfortran

MEXA64Files = $(O)/BKMP2Mex.mexa64 $(O)/H3LSTHMex.mexa64

OBJS = $(O)/bkmp2.o $(O)/lsth_h3.o

#QMLibs = $(O)/libqmdyn.a

.DEFAULT_GOAL := $(O)/H3LSTHMex.mexa64

all: $(MEXA64Files)

#$(EXE) : $(OBJS)
#	$(Link) $(CXXFLAGS) -o $(EXE) $(OBJS) $(LIBS)

$(O)/%.o: %.c
	cd $(O) ; $(CC)  $(cFLAGS) -c $(SRC)/$<
$(O)/%.o: %.C
	cd $(O) ; $(CXX) $(CXXFLAGS) -c $(SRC)/$<
$(O)/%.o: %.F
	cd $(O) ; $(F77) $(FFLAGS) -c $(SRC)/$<
$(O)/%.o: %.F90
	cd $(O) ; $(F90) $(FFLAGS) -c $(SRC)/$<
%io.C: %.h
	perl io.pl $<

$(QMLibs): $(OBJS)
	cd $(O); ar -crusv $(QMLibs) $(OBJS)

$(O)/%.mexa64: $(O)/%.o $(OBJS) #$(QMLibs)
	cd $(O); $(Link) -shared $(CXXFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o *~ *.mod $(EXE) depend $(MEXA64Files) $(QMLibs) $(OBJS)

depend :
	$(CXX) $(CXXFLAGS) -MM *.[cC] | perl dep.pl | tee $@
	sh fortran-dep.sh *.F >> $@	

ifneq ($(MAKECMDGOALS), clean)
include depend
endif
