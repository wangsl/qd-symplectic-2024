
# $Id: Makefile,v 1.1 2005/07/18 20:32:28 wangsl Exp wangsl $

.DEFAULT: .F .For .c .C .cpp .cc .F90
.SUFFIXES: .F .For .c .C .cpp .cc .F90

SRC = $(PWD)
O = $(SRC)/../../build/common

F77 = gfortran
CC = gcc
CCC = g++

CFLAGS = -O3 -fPIC -fopenmp -I$(MATLAB_ROOT)/extern/include

FFLAGS = -O3 -fPIC -fopenmp -fallow-argument-mismatch

Link = $(CCC) $(FFLAGS)

LIBS = #-lifcoremt -lifport

MEXA64Files = $(O)/GaussLegendreGridsMex.mexa64 $(O)/fcnpakMex.mexa64

#.DEFAULT_GOAL := $(MEXA64Files)

all: $(MEXA64Files)

OBJS = $(O)/gaussq.o  $(O)/fcnpak.o  $(O)/dgamma.o

#EXENAME =

#$(EXENAME) : $(OBJS) 
#	$(Link) -o $(EXENAME) $(OBJS) $(LIBS)\

$(O)/%.mexa64: $(O)/%.o $(OBJS) 
	cd $(O); $(Link) -shared $(CXXFLAGS) -o $@ $^ $(LIBS)

$(O)/%.o: %.c
	cd $(O) ; $(CC) $(CFLAGS) -c $(SRC)/$<
$(O)/%.o: %.cc
	cd $(O) ; $(CCC) $(CFLAGS) -c $(SRC)/$<
$(O)/%.o: %.cpp
	cd $(O) ; $(CCC) $(CFLAGS) -c $(SRC)/$<
$(O)/%.o: %.C
	cd $(O) ; $(CCC) $(CFLAGS) -c $(SRC)/$<
$(O)/%.o: %.F
	cd $(O) ; $(F77) $(FFLAGS) -c $(SRC)/$<
$(O)/%.o: %.for
	cd $(O) ; $(F77) $(FFLAGS) -c $(SRC)/$<
$(O)/%.o: %.f90
	cd $(O) ; $(F90) $(FFLAGS) -c $(SRC)/$<

clean:
	rm -f $(O)/*.o *~ $(EXENAME) $(O)/*.mod *io.C depend $(O)/*.mexa64 $(O)/*.optrpt

.PRECIOUS: %io.C %.C %.c %.F %.h $(O)/%.o

depend :
	$(CCC) $(CFLAGS) -MM *.[cC] | perl dep.pl > $@
	bash dep.sh *.F >> $@

ifneq ($(MAKECMDGOALS), clean)
include depend
endif

