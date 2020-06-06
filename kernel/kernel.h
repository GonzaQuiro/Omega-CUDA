/*
 * kernel.h
 *
 *  Created on: 01/07/2019
 *      Author: Gonzalo Damian Quiroga
 */

#ifndef KERNEL_H_
#define KERNEL_H_

__global__ void rk4(double m1, double m2, double *tau, double *x0, double *x1, double *x2, double *x3, double *x4, double *x5, double *x6, double *x7, double *x8, double *x9, double *x10, int n_ic, double h, int MainLoop);
__global__ void dp45(double m1, double m2, double *tau, double *x0, double *x1, double *x2, double *x3, double *x4, double *x5, double *x6, double *x7, double *x8, double *x9, double *x10, int n_ic, double h, double tol, int MainLoop);
__global__ void rkf78(double m1, double m2, double *tau, double *x0, double *x1, double *x2, double *x3, double *x4, double *x5, double *x6, double *x7, double *x8, double *x9, double *x10, int n_ic, double h, double tol, int MainLoop);

#endif /* KERNEL_H_ */
