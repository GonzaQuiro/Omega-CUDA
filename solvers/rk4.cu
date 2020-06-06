// Version 1.0.0 CUDA-C: Runge-Kutta 4 method
// Dr. Gonzalo Damián Quiroga
#include "../ODEs/rhs.h"

void rk4(double m1, double m2, double *tau, double *x0, double *x1, double *x2, double *x3, double *x4, double *x5, double *x6, double *x7, double *x8, double *x9, double *x10, int n_ic, double h, int MainLoop)
{
	const int id = blockIdx.x*blockDim.x + threadIdx.x;
	double k1[11], k2[11], k3[11], k4[11];
	double dx[11];
	int count;

	if (id < n_ic)
	{	
		//Stop Condition
		if(x0[id]<0.05)
		{
			//Coefficinet k1:
			rhs(m1, m2, x0[id], x1[id], x2[id], x3[id], x4[id], x5[id], x6[id], x7[id], x8[id], x9[id], x10[id], 
			&dx[0], &dx[1], &dx[2], &dx[3], &dx[4], &dx[5], &dx[6], &dx[7], &dx[8], &dx[9], &dx[10]);
			
			for(count=0;count<11;count++)
				k1[count] = h*dx[count];
				
			//Coefficinet k2:
			rhs(m1, m2, x0[id] + 0.5*k1[0], x1[id] + 0.5*k1[1], x2[id] + 0.5*k1[2], 
			x3[id] + 0.5*k1[3],x4[id] + 0.5*k1[4], x5[id] + 0.5*k1[5], x6[id] + 0.5*k1[6], x7[id] + 0.5*k1[7],
			x8[id] + 0.5*k1[8],x9[id] + 0.5*k1[9], x10[id] + 0.5*k1[10],
			&dx[0], &dx[1], &dx[2], &dx[3], &dx[4], &dx[5], &dx[6], &dx[7], &dx[8], &dx[9], &dx[10]);
			
			for(count=0;count<11;count++)
				k2[count] = h*dx[count];
			
			//Coefficinet k3:
			rhs(m1, m2, x0[id] + 0.5*k2[0], x1[id] +0.5*k2[1], x2[id]+0.5*k2[2], 
			x3[id]+0.5*k2[3], x4[id]+0.5*k2[4], x5[id]+0.5*k2[5], x6[id]+0.5*k2[6], x7[id]+0.5*k2[7],
			x8[id]+0.5*k2[8], x9[id]+0.5*k2[9], x10[id]+0.5*k2[10],  
			&dx[0], &dx[1], &dx[2], &dx[3], &dx[4], &dx[5], &dx[6], &dx[7], &dx[8], &dx[9], &dx[10]);

			for(count=0;count<11;count++)
				k3[count] = h*dx[count];

			//Coefficinet k4:
			rhs(m1, m2, x0[id]+k3[0],x1[id]+k3[1], x2[id]+k3[2], x3[id]+k3[3],
			x4[id]+k3[4], x5[id]+k3[5],x6[id]+k3[6], x7[id]+k3[7], 
			x8[id]+k3[8],x9[id]+k3[9], x10[id]+k3[10], 
			&dx[0], &dx[1], &dx[2], &dx[3], &dx[4], &dx[5], &dx[6], &dx[7], &dx[8], &dx[9], &dx[10]);

			for(count=0;count<11;count++)
				k4[count] = h*dx[count];

			tau[id]=MainLoop*h;
			x0[id]=x0[id]+(k1[0]+2.0*k2[0]+2.0*k3[0]+k4[0])/6.0;
			x1[id]=x1[id]+(k1[1]+2.0*k2[1]+2.0*k3[1]+k4[1])/6.0;
			x2[id]=x2[id]+(k1[2]+2.0*k2[2]+2.0*k3[2]+k4[2])/6.0;
			x3[id]=x3[id]+(k1[3]+2.0*k2[3]+2.0*k3[3]+k4[3])/6.0;
			x4[id]=x4[id]+(k1[4]+2.0*k2[4]+2.0*k3[4]+k4[4])/6.0;
			x5[id]=x5[id]+(k1[5]+2.0*k2[5]+2.0*k3[5]+k4[5])/6.0;
			x6[id]=x6[id]+(k1[6]+2.0*k2[6]+2.0*k3[6]+k4[6])/6.0;
			x7[id]=x7[id]+(k1[7]+2.0*k2[7]+2.0*k3[7]+k4[7])/6.0;
			x8[id]=x8[id]+(k1[8]+2.0*k2[8]+2.0*k3[8]+k4[8])/6.0;
			x9[id]=x9[id]+(k1[9]+2.0*k2[9]+2.0*k3[9]+k4[9])/6.0;
			x10[id]=x10[id]+(k1[10]+2.0*k2[10]+2.0*k3[10]+k4[10])/6.0;
		}
		else
		{
			return;
		}
	}
}

