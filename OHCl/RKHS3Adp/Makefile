
# $Id$

.DEFAULT: .F .F90 .c .C
.SUFFIXES: .F .F90 .c .C

SRC = $(PWD)
base_dir=$(SRC)/..
O = $(base_dir)/build/OHCl/RKHS3Adp

F77 = ifort
F90 = ifort
CC = icc
CXX = icpc

CFLAGS =  

CXXFLAGS = -std=c++11 -I$(MATLAB_ROOT)/extern/include -I$(base_dir)/common-cuda

FFLAGS = -132 -I$(MATLAB_ROOT)/extern/include

Link = $(CXX)

LIBS = -lifcoremt -L$(base_dir)/build/common-cuda -lqmdyn

MEXA64Files = $(O)/OHClRKHSMex.mexa64

OBJS = $(O)/OHClRKHSMex.o $(O)/rkhs_pes_3Adp.o

.DEFAULT_GOAL := $(O)/OHClRKHSMex.mexa64

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

$(O)/%.mexa64: $(O)/%.o $(OBJS)
	cd $(O); $(Link) -shared $(CXXFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o *~ *.mod $(EXE) depend $(MEXA64Files) $(QMLibs) $(OBJS)

depend :
	$(CXX) $(CXXFLAGS) -MM *.[cC] | perl dep.pl | tee $@
	sh fortran-dep.sh *.F >> $@

ifneq ($(MAKECMDGOALS), clean)
include depend
endif

