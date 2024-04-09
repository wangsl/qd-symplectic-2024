#!/bin/bash

shopt -s expand_aliases
alias die='_error "Error in file $0 at line $LINENO:"'
alias warn='_warn "Warning in file $0 at line $LINENO:"'

function setup_cuda_compilers()
{
    local version=11.0.2
    module load cuda/$version
    export NVCC_PATH=$INTEL_WRAPPER_PATH/cuda/${version}/bin
    export NVCC_BIN_PATH=$(dirname $(which nvcc))
}

function special_rules()
{
    return
    
    local arg=
    for arg in "$@"; do
	echo $arg
    done
}

function main()
{
    source /opt/apps/lmod/lmod/init/bash
    export LMOD_DISABLE_SAME_NAME_AUTOSWAP=yes
    module unuse /share/apps/modulefiles
    module use /ext3/apps/modulefiles   
    module purge

    export PKG_CONFIG_PATH=
    export CPATH=
    export LIBRARY_PATH=
    export LD_LIBRARY_PATH=

    module load intel/19.1.1
    module load matlab/2020a

    export INTEL_WRAPPER_PATH="/ext3/apps/utils/intel"

    setup_cuda_compilers

    local util=$INTEL_WRAPPER_PATH/util.bash
    if [ -e $util ]; then source $util; fi
    
    export SPECIAL_RULES_FUNCTION=special_rules
    if [ "$SPECIAL_RULES_FUNCTION" != "" ]; then
	export BUILD_WRAPPER_SCRIPT=$(readlink -e $0)
    fi

    export GNU_BIN_PATH=$(dirname $(which gcc))
    export INTEL_BIN_PATH=$(dirname $(which icc))
    #export INTEL_MPI_BIN_PATH=$(dirname $(which mpicc))

    export INVALID_FLAGS="-O -O0 -O1 -O2 -O3 -g -g0"
    
    export INVALID_FLAGS_FOR_GNU_COMPILERS=""
    export OPTIMIZATION_FLAGS_FOR_GNU_COMPILERS="-fPIC -fopenmp -mavx2"
    
    export INVALID_FLAGS_FOR_INTEL_COMPILERS="-lm -xhost -fast"
    
    export OPTIMIZATION_FLAGS_FOR_INTEL_COMPILERS="-fPIC -unroll -ip -axCORE-AVX512 -qopenmp -qopt-report-stdout -qopt-report-phase=openmp"
    
    export OPTIMIZATION_FLAGS_FOR_INTEL_FORTRAN_COMPILERS="-fPIC -unroll -ip -axCORE-AVX512 -qopenmp -qopt-report-stdout -qopt-report-phase=openmp"

    #export INVALID_FLAGS_FOR_NVCC_COMPILERS=""
    #export OPTIMIZATION_FLAGS_FOR_NVCC_COMPILERS="-Wno-deprecated-gpu-targets"
    
    export OPTIMIZATION_FLAGS="-O3"
    
    export CPPFLAGS=$(for inc in $(env -u INTEL_INC -u MKL_INC | grep _INC= | cut -d= -f2); do echo '-I'$inc; done | xargs)
    export LDFLAGS=$(for lib in $(env | grep _LIB= | cut -d= -f2); do echo '-L'$lib; done | xargs)
    
    prepend_to_env_variable INCLUDE_FLAGS "$CPPFLAGS"
    prepend_to_env_variable LINK_FLAGS "$LDFLAGS"
    
    export INCLUDE_FLAGS_FOR_INTEL_COMPILERS="-I$INTEL_INC -I$MKL_INC"
    
    export LINK_FLAGS_FOR_INTEL_COMPILERS="-shared-intel"
    export EXTRA_LINK_FLAGS="$(LD_LIBRARY_PATH_to_rpath)"
    
    if [ "$DEBUG_LOG_FILE" != "" ]; then
	echo "" | tee $DEBUG_LOG_FILE
    fi
    
    export LD_RUN_PATH=$LD_LIBRARY_PATH

    local prefix=$(pwd)
    if [ "$prefix" == "" ]; then
	local dir=$(readlink -e $(dirname $0))
	dir="$dir/local"
	if [ -d $dir ]; then prefix=$dir; fi
    fi
    if [ "$prefix" == "" ]; then die "no prefix defined"; fi

    #export N_MAKE_THREADS=60
    #export DEFAULT_COMPILER="GNU"
    
    local args="$@"
    local arg=
    for arg in $args; do
	
	case $arg in
	    
	    configure|conf)
		echo " Run configuration ..."
		export PATH=.:$INTEL_WRAPPER_PATH:$NVCC_PATH:$PATH
		
		if [ "$DEFAULT_COMPILER" != "GNU" ]; then
		    export CC=icc
                    export CXX=icpc
                    export FC=ifort
		    export F77=ifort
		fi
		
		./configure --build=x86_64-centos-linux \
			    --prefix=$prefix
		;;
	    
	    cmake)
		module load cmake/gcc/3.17.3
		export PATH=.:$INTEL_WRAPPER_PATH:$PATH
		
		export CMAKE_INCLUDE_PATH=$(env | grep _INC= | cut -d= -f2 | xargs | sed -e 's/ /:/g')
		export CMAKE_LIBRARY_PATH=$(env | grep _LIB= | cut -d= -f2 | xargs | sed -e 's/ /:/g')
		
                export CC=icc
                export CXX=icpc
		cmake \
		    -DCMAKE_BUILD_TYPE=release \
                    -DBUILD_SHARED_LIBS::BOOL=ON \
                    -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
                    -DCMAKE_SKIP_RPATH:BOOL=ON \
		    -DCMAKE_INSTALL_PREFIX:PATH=$prefix \
		    ../breakdancer
                ;;
	    
	    make)
		export PATH=.:$INTEL_WRAPPER_PATH:$NVCC_PATH:$PATH
		echo " Run make"
		eval "$args" 
		exit
		;;

	    a2so)
		export PATH=.:$INTEL_WRAPPER_PATH:$PATH
		icc -shared -o libsuitesparse.so  \
		    -Wl,--whole-archive \
		    libamd.a \
		    -Wl,--no-whole-archive \
		    -L$MKL_ROOT/lib/intel64 -lmkl_rt
		exit
		;;
	    
	    *)
		die "Usage: $0 <argument>: configure make"
		;;
	esac

	args=$(eval "echo $args | sed -e 's/$arg //'")
    done
}

#############################################################
# do the main work here, do not modify the follwoing part,  #
# we just need to modify the main function                  #
#############################################################

if [ "$TO_SOURCE_BUILD_WRAPPER_SCRIPT" == "" ]; then
    main "$@"
    exit
else
    unset -f main
fi

#############################################################
# End here, do not add anything after this line             #
#############################################################
