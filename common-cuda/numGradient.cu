
#include "cudaUtils.h"
#include "numGradient.h"
#include "complex.h"

__constant__ double gradient_coeffients_dev[_GradientMaxSize_];

void Num1ststGradient::copy_gradient_coefficients_to_device(const int n_points)
{
  const int n = n_points/2;
  insist(n_points == 2*n+1);

  size_t size = 0;
  checkCudaErrors(cudaGetSymbolSize(&size, gradient_coeffients_dev));
  insist(size > n*sizeof(double));
  
  const double *coefs = coefficients(n_points);
  
  checkCudaErrors(cudaMemcpyToSymbolAsync(gradient_coeffients_dev, coefs, n*sizeof(double),
					  0, cudaMemcpyHostToDevice));
}

template<class T> __global__ void gradients_1_3d(const int n1, const int n2, const int n3,
						 const int n1_surface, const int n_points,
						 const double dx1,
						 const T *f, T *vals, T *grads)
{
  const int index = threadIdx.x + blockDim.x*blockIdx.x;
  
  const int n = n_points/2;
  assert(n_points == 2*n+1);
  
  if(index < n2*n3) {
    int j = -1; int k = -1;
    cudaUtils::index_2_ij(index, n2, n3, j, k);
    
    vals[index] = f[cudaUtils::ijk_2_index(n1, n2, n3, n1_surface, j, k)];
    
    T g = 0.0;
    for(int ii = 0; ii < n; ii++) {
      g += gradient_coeffients_dev[ii]*
	(f[cudaUtils::ijk_2_index(n1, n2, n3, n1_surface+ii+1, j, k)] - 
	 f[cudaUtils::ijk_2_index(n1, n2, n3, n1_surface-ii-1, j, k)]);
    }
    grads[index] = g/dx1;
  }
}

template<class T> __global__ void gradients_2_3d(const int n1, const int n2, const int n3,
						 const int n2_surface, const int n_points,
						 const double dx2,
						 const T *f, T *vals, T *grads)
{
  const int index = threadIdx.x + blockDim.x*blockIdx.x;

  const int n = n_points/2;
  assert(n_points == 2*n+1);
  
  if(index < n1*n3) {
    int i = -1; int k = -1;
    cudaUtils::index_2_ij(index, n1, n3, i, k);
    
    vals[index] = f[cudaUtils::ijk_2_index(n1, n2, n3, i, n2_surface, k)];
    
    T g = 0.0;
    for(int ii = 0; ii < n; ii++) {
      g += gradient_coeffients_dev[ii]*
	(f[cudaUtils::ijk_2_index(n1, n2, n3, i, n2_surface+ii+1, k)] - 
	 f[cudaUtils::ijk_2_index(n1, n2, n3, i, n2_surface-ii-1, k)]);
    }
    grads[index] = g/dx2;
  }
}

template __global__ void gradients_1_3d<double>(const int n1, const int n2, const int n3,
						const int n1_surface, const int n_points,
						const double dx1,
						const double *f, double *vals, double *grads);

template __global__ void gradients_1_3d<Complex>(const int n1, const int n2, const int n3,
						 const int n1_surface, const int n_points,
						 const double dx1,
						 const Complex *f, Complex *vals, Complex *grads);

template __global__ void gradients_2_3d<double>(const int n1, const int n2, const int n3,
						const int n2_surface, const int n_points,
						const double dx2,
						const double *f, double *vals, double *grads);

template __global__ void gradients_2_3d<Complex>(const int n1, const int n2, const int n3,
						 const int n2_surface, const int n_points,
						 const double dx2,
						 const Complex *f, Complex *vals, Complex *grads);

void Num1ststGradient::gradients_1_3d(const int n1, const int n2, const int n3,
				      const int n1_surface, const int n_points,
				      const double dx1,
				      const double *f, double *values, double *grads)
{
  const int n_threads = _NTHREADS_;
  int n_blocks = cudaUtils::number_of_blocks(n_threads, n2*n3);
  
  ::gradients_2_3d<double><<<n_blocks, n_threads>>>
      (n1, n2, n3, n1_surface, n_points, dx1, f, values, grads);
}

void Num1ststGradient::gradients_2_3d(const int n1, const int n2, const int n3,
				      const int n2_surface, const int n_points,
				      const double dx2,
				      const double *f, double *values, double *grads)
{
  const int n_threads = _NTHREADS_;
  int n_blocks = cudaUtils::number_of_blocks(n_threads, n1*n3);
  
  ::gradients_2_3d<double><<<n_blocks, n_threads>>>
      (n1, n2, n3, n2_surface, n_points, dx2, f, values, grads);
}
