
#ifndef CUDA_UTILS_H
#define CUDA_UTILS_H

#include <cstdlib>
#include <cstdio>

#include <chrono>
#include <thread>
#include <hipblas.h>
#include <hipfft.h>

#define CHECK(cmd)                              \
  {                                             \
    hipError_t error  = cmd;                   \
    if (error != hipSuccess) {                                         \
      fprintf(stderr, "error: '%s'(%d) at %s:%d\n", hipGetErrorString(error), error,__FILE__, __LINE__); \
      exit(EXIT_FAILURE);                                               \
    }                                                                   \
  }

#define checkCudaErrors(x) CHECK(x)

#define _CUDA_FREE_(x) if(x) { checkCudaErrors(hipFree(x)); x = 0; }

#define _CUDA_STREAM_CREATE_(x) {					\
    if(!x) {								\
      x = (hipStream_t *) malloc(sizeof(hipStream_t));		\
      insist(x);							\
      checkCudaErrors(hipStreamCreate(x));				\
    }									\
  }			       

#define _CUDA_EVENT_CREATE_(x) {					\
    if(!x) {								\
      x = (hipEvent_t *) malloc(sizeof(hipEvent_t));			\
      insist(x);							\
      checkCudaErrors(hipEventCreateWithFlags(x, hipEventDisableTiming)); \
    }									\
  }

#define _CUDA_STREAM_DESTROY_(x) {		\
    if(x) {					\
      checkCudaErrors(hipStreamDestroy(*x));	\
      free(x);					\
      x = 0;					\
    }						\
  }

#define _CUDA_EVENT_DESTROY_(x) {		\
    if(x) {					\
      checkCudaErrors(hipEventDestroy(*x));	\
      free(x);					\
      x = 0;					\
    }						\
  }

#define _NTHREADS_ 512
#define _POTENTIAL_CUTOFF_ -1.0e+6

namespace cudaUtils {

  __device__ __host__ inline int ij_2_index(const int n1, const int n2, const int i, const int j)
  { return j*n1 + i; }
  
  __device__ __host__ inline void index_2_ij(const int index, const int n1, const int n2, int &i, int &j)
  {  j = index/n1; i = index - j*n1; }
  
  __device__ __host__ inline int ijk_2_index(const int n1, const int n2, const int n3, 
					     const int i, const int j, const int k)
  { return (k*n2 + j)*n1 + i; }
  
  __device__ __host__ inline void index_2_ijk(const int index, const int n1, const int n2, const int n3, 
					      int &i, int &j, int &k)
  {
    int ij = -1;
    index_2_ij(index, n1*n2, n3, ij, k);
    index_2_ij(ij, n1, n2, i, j);
  }
  
  inline int number_of_blocks(const int n_threads, const int n)
  { return n/n_threads*n_threads == n ? n/n_threads : n/n_threads+1; }

  inline int n_devices()
  {
    int _n_devices = -1;
    checkCudaErrors(hipGetDeviceCount(&_n_devices));
    return _n_devices;
  }
  
  void device_memory_usage();

  void cufft_work_size(const hipfftHandle &plan, const char *type = 0);
}

inline char *time_now()
{
  std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
  std::time_t time_now  = std::chrono::system_clock::to_time_t(now);
  char *time = std::ctime(&time_now); 
  char *pch = strchr(time, '\n');
  if(pch) pch[0] = '\0';
  return time;
}

#endif /* CUDA_UTILS_H */
  

