
#ifndef EVOLUTION_UTILS_H
#define EVOLUTION_UTILS_H

#if defined __HIPCC__
#include <hip/hip_runtime.h>
#endif

#include "matlabStructures.h"
#include "numGradient.h"

#define _RealPart_ 1
#define _ImagPart_ 2

#define _RotStatesAll_ 0
#define _RotStatesOdd_ 1
#define _RotStatesEven_ 2

#define _DumpMaxSize_ 2048
#define _EnergiesMaxSize_ 1024

namespace EvolutionUtils {

#if defined __NVCC__ || defined __HIPCC__
  
  struct RadialCoordinate
  {
    double left;
    double dr;
    double mass;
    double dump[_DumpMaxSize_];
    int n;
  };
  
  inline void copy_radial_coordinate_to_device(const RadialCoordinate &r_dev, 
					       const ::RadialCoordinate *r)
  {
    size_t size = 0;
    checkCudaErrors(cudaGetSymbolSize(&size, r_dev));
    insist(size == sizeof(RadialCoordinate));

    RadialCoordinate r_;
    r_.left = r->left;
    r_.dr = r->dr;
    r_.mass = r->mass;
    r_.n = r->n;
    
    insist(r->n <= _DumpMaxSize_);

    memset(r_.dump, 0, _DumpMaxSize_*sizeof(double));
    memcpy(r_.dump, r->dump, r->n*sizeof(double));
    
    checkCudaErrors(cudaMemcpyToSymbolAsync(r_dev, &r_, 
					    sizeof(RadialCoordinate), 0, 
					    cudaMemcpyHostToDevice));
  }
  
#endif
  
  inline int int_to_even_right(const int i) { return i/2*2 == i ? i : i+1; } 
  inline int int_to_even_left(const int i) { return i/2*2 == i ? i : i-1; }  
  inline int int_to_odd_right(const int i) { return i/2*2 == i ? i+1 : i; }
  inline int int_to_odd_left(const int i) { return i/2*2 == i ? i-1 : i; }

  inline double au_to_fs(const double t) { return t*0.0241888425; }
  inline double fs_to_au(const double t) { return t/au_to_fs(1.0); }
}

#if defined __NVCC__ || defined __HIPCC__
extern __constant__ EvolutionUtils::RadialCoordinate r1_dev;
extern __constant__ EvolutionUtils::RadialCoordinate r2_dev;
extern __constant__ double energies_dev[_EnergiesMaxSize_];
extern __constant__ double gradient_coeffients_dev[_GradientMaxSize_];
extern __constant__ double potential_cutoff;
//extern __constant__ double kinetic_cutoff;
#endif

#endif /* EVOLUTION_UTILS_H */


