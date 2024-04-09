
# $Id$

.DEFAULT: .F .F90 .c .C
.SUFFIXES: .F .F90 .c .C

SRC = $(PWD)
base_dir = $(SRC)/..
O = $(base_dir)/build/common-cuda

F77 = ifort
F90 = ifort
CC = icc
CXX = icpc
NVCC = nvcc

CFLAGS =

FFLAGS = -132

NVCCINCLUDE = -I$(CUDA_ROOT)/samples/common/inc \
	-I$(MATLAB_ROOT)/extern/include

CUDAArchs= \
	-gencode arch=compute_20,code=sm_20 \
        -gencode arch=compute_30,code=sm_30 \
        -gencode arch=compute_35,code=sm_35 \
        -gencode arch=compute_37,code=sm_37 \
        -gencode arch=compute_50,code=sm_50 \
	-gencode arch=compute_52,code=sm_52 \
        -gencode arch=compute_52,code=compute_52 \
	-Xcompiler=\"-fPIC -pthread -fexceptions -m64 -fopenmp\"

CUDAArchs= \
	-gencode arch=compute_37,code=sm_37 \
	-gencode arch=compute_52,code=sm_52 \
	-gencode arch=compute_52,code=compute_52 \
	-Xcompiler=\"-fPIC -pthread -fexceptions -m64 -fopenmp\"

CUDAArchs= \
	-gencode arch=compute_37,code=sm_37 \
	-gencode arch=compute_52,code=sm_52 \
	-gencode arch=compute_52,code=compute_52 \
	-Xcompiler=\"-fPIC -pthread -fexceptions -m64 -fopenmp\"


NVCCFLAGS = $(NVCCINCLUDE) $(CUDAArchs) -prec-div=true -prec-sqrt=true -rdc=true -std=c++11 

CXXFLAGS = -std=c++11 $(NVCCINCLUDE)

Link = $(CXX)

LIBS = -lifcoremt -L$(CUDA_LIB) -lcufft -lcublas -lcudart #-lmkl_rt

MEXA64Files = $(O)/cudaSymplectic.mexa64 #$(O)/DMBEIVMex.mexa64

CUDAObjs = $(O)/cudaOpenmpMD.o $(O)/cudaUtils.o $(O)/wavepacketson1device.o \
	$(O)/omegawavepacket.o $(O)/numGradient.o  $(O)/reactProb.o
 
CUDALinkObj = $(O)/cudalink.o

OBJS = $(O)/cudaSymplectic.o  $(O)/matlabUtils.o \
	$(O)/matlabStructures.o  $(O)/matlabStructuresio.o \
	$(O)/matlabData.o \
	$(O)/die.o  $(O)/indent.o  $(O)/out.o \
	$(O)/rmatalgo.o  $(O)/rmat.o  $(O)/rmato.o \
	$(O)/mex2matlab.o \
	$(CUDAObjs) $(CUDALinkObj)

QMLibs = $(O)/libqmdyn.a

#.DEFAULT_GOAL := $(O)/DMBEIVMex.mexa64
.DEFAULT_GOAL := $(O)/cudaSymplectic.mexa64

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
$(O)/%.o: %.cu
	cd $(O) ; $(NVCC) $(NVCCFLAGS) -dc $(SRC)/$<

$(CUDALinkObj): $(CUDAObjs)
	cd $(O); $(NVCC) $(CUDAArchs) -dlink $(CUDAObjs) -o $(CUDALinkObj)

%io.C: %.h
	perl io.pl $<

$(QMLibs): $(OBJS)
	cd $(O); ar -crusv $(QMLibs) $(OBJS)

$(O)/%.mexa64: $(O)/%.o $(QMLibs)
	cd $(O); $(Link) -shared $(CXXFLAGS) -o $@ $^ $(LIBS)

clean:
	rm -f *.o *~ *.mod $(EXE) depend $(MEXA64Files) $(QMLibs) $(OBJS) $(O)/*.optrpt

cleancuda:
	rm -rf $(CUDAObjs) $(CUDALinkObj)

depend :
	$(CXX) $(CXXFLAGS) -MM *.[cC] | perl dep.pl > $@
	$(NVCC) $(CXXFLAGS) -M *.cu | perl dep.pl >> $@
	sh dep.sh *.F >> $@

ifneq ($(MAKECMDGOALS), clean)
include depend
endif

