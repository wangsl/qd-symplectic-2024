
#include <iostream>
#include <iomanip>

#include <hip/hip_runtime.h>
#include <hipfft.h>

#include "cudaUtils.h"
#include "matlabUtils.h"

void cudaUtils::device_memory_usage()
{
  std::ios_base::fmtflags old_flags = std::cout.flags();
  std::streamsize old_precision = std::cout.precision();
  
  std::cout.precision(2);
  
  int device_index = -1;
  checkCudaErrors(hipGetDevice(&device_index));
  
  size_t free_byte = 0;
  size_t total_byte = 0;
  checkCudaErrors(hipMemGetInfo(&free_byte, &total_byte));
  
  std::cout << " Device: " << device_index
	    << " GPU memory usage:" 
	    << " used = " << std::fixed << (total_byte-free_byte)/1024.0/1024.0 << "MB,"
	    << " free = " << free_byte/1024.0/1024.0 << "MB,"
	    << " total = " << total_byte/1024.0/1024.0 << "MB" << std::endl;

  std::cout.flags(old_flags);
  std::cout.precision(old_precision);
}

void cudaUtils::cufft_work_size(const hipfftHandle &plan, const char *type)
{
  std::ios_base::fmtflags old_flags = std::cout.flags();
  std::streamsize old_precision = std::cout.precision();
  
  std::cout.precision(2);
  
  size_t cufft_work_size = 0;
  insist(hipfftGetSize(plan, &cufft_work_size) == HIPFFT_SUCCESS);
  std::cout << " Work area size for CUFFT";
  if(type) std::cout << " " << type;
  std::cout << ":" << " " << std::fixed << cufft_work_size/1024.0/1024.0 << "MB" << std::endl;

  std::cout.flags(old_flags);
  std::cout.precision(old_precision);
}

