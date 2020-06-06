// Version 1.0.0 CUDA-C: Omega
// Dr. Gonzalo Damián Quiroga
// Universidad Nacional de Córdoba
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include "../ODEs/rhs.h"

__global__ 
#include "../solvers/rk4.cu"

__global__ 
#include "../solvers/dp45.cu"

__global__ 
#include "../solvers/rkf78.cu"

__device__ 
#include "../ODEs/rhs.c"

