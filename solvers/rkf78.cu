// Version 1.0.0 CUDA-C: Runge-Kutta 4 method
// Dr. Gonzalo Damián Quiroga
#include "../ODEs/rhs.h"

void rkf78(double m1, double m2, double *tau, double *x0, double *x1, double *x2, double *x3, double *x4, double *x5, double *x6, double *x7, double *x8, double *x9, double *x10, int n_ic, double h, double tol, int MainLoop)
{
	const int id = blockIdx.x*blockDim.x + threadIdx.x;
	double k1[11], k2[11], k3[11], k4[11], k5[11], k6[11], k7[11], k8[11], k9[11], k10[11], k11[11], k12[11], k13[11];
	double dx[11], delta[11];
	double error;
	int count;
	int end=1;
	int loop=0;

	if (id < n_ic)
	{	
		//Stop Condition
		if(x0[id]<0.1)
		{
			while(loop<end)
			{
				start:

				//Coefficinet k1:
				rhs(m1, m2, x0[id], x1[id], x2[id], x3[id], x4[id], x5[id], x6[id], x7[id], x8[id], x9[id], x10[id], 
				&dx[0], &dx[1], &dx[2], &dx[3], &dx[4], &dx[5], &dx[6], &dx[7], &dx[8], &dx[9], &dx[10]);
				
				for(count=0;count<11;count++)
					k1[count] = h*dx[count];
					
				//Coefficinet k2:
				for(count=0;count<11;count++)
					delta[count] = 2.0*k1[count]/27.0;

				rhs(m1, m2, 
				x0[id] + delta[0], 
				x1[id] + delta[1],
				x2[id] + delta[2], 
				x3[id] + delta[3], 
				x4[id] + delta[4], 
				x5[id] + delta[5], 
				x6[id] + delta[6], 
				x7[id] + delta[7],
				x8[id] + delta[8], 
				x9[id] + delta[9], 
				x10[id] + delta[10],
				&dx[0], &dx[1], &dx[2], &dx[3], &dx[4], &dx[5], &dx[6], &dx[7], &dx[8], &dx[9], &dx[10]);
				
				for(count=0;count<11;count++)
					k2[count] = h*dx[count];
				
				//Coefficinet k3:
				for(count=0;count<11;count++)
					delta[count] = k1[count]/36.0 + k2[count]/12.0;

				rhs(m1, m2, 
				x0[id] + delta[0], 
				x1[id] + delta[1],
				x2[id] + delta[2], 
				x3[id] + delta[3], 
				x4[id] + delta[4], 
				x5[id] + delta[5], 
				x6[id] + delta[6], 
				x7[id] + delta[7],
				x8[id] + delta[8], 
				x9[id] + delta[9], 
				x10[id] + delta[10],
				&dx[0], &dx[1], &dx[2], &dx[3], &dx[4], &dx[5], &dx[6], &dx[7], &dx[8], &dx[9], &dx[10]);

				for(count=0;count<11;count++)
					k3[count] = h*dx[count];
				
				//Coefficinet k4:
				for(count=0;count<11;count++)
					delta[count] = k1[count]/24.0+k3[count]/8.0;

				rhs(m1, m2, 
				x0[id] + delta[0], 
				x1[id] + delta[1],
				x2[id] + delta[2], 
				x3[id] + delta[3], 
				x4[id] + delta[4], 
				x5[id] + delta[5], 
				x6[id] + delta[6], 
				x7[id] + delta[7],
				x8[id] + delta[8], 
				x9[id] + delta[9], 
				x10[id] + delta[10],
				&dx[0], &dx[1], &dx[2], &dx[3], &dx[4], &dx[5], &dx[6], &dx[7], &dx[8], &dx[9], &dx[10]);

				for(count=0;count<11;count++)
					k4[count] = h*dx[count];

				//Coefficinet k5:
				for(count=0;count<11;count++)
					delta[count] = 5.0*k1[count]/12.0-25.0*k3[count]/16.0+25.0*k4[count]/16.0;

				rhs(m1, m2, 
				x0[id] + delta[0], 
				x1[id] + delta[1],
				x2[id] + delta[2], 
				x3[id] + delta[3], 
				x4[id] + delta[4], 
				x5[id] + delta[5], 
				x6[id] + delta[6], 
				x7[id] + delta[7],
				x8[id] + delta[8], 
				x9[id] + delta[9], 
				x10[id] + delta[10],
				&dx[0], &dx[1], &dx[2], &dx[3], &dx[4], &dx[5], &dx[6], &dx[7], &dx[8], &dx[9], &dx[10]);
			
				for(count=0;count<11;count++)
					k5[count] = h*dx[count];

				//Coefficinet k6:
				for(count=0;count<11;count++)
					delta[count] = k1[count]/20.0+k4[count]/4.0+k5[count]/5.0;

				rhs(m1, m2, 
				x0[id] + delta[0], 
				x1[id] + delta[1],
				x2[id] + delta[2], 
				x3[id] + delta[3], 
				x4[id] + delta[4], 
				x5[id] + delta[5], 
				x6[id] + delta[6], 
				x7[id] + delta[7],
				x8[id] + delta[8], 
				x9[id] + delta[9], 
				x10[id] + delta[10],
				&dx[0], &dx[1], &dx[2], &dx[3], &dx[4], &dx[5], &dx[6], &dx[7], &dx[8], &dx[9], &dx[10]);
			
				for(count=0;count<11;count++)
					k6[count] = h*dx[count];

				//Coefficinet k7:
				for(count=0;count<11;count++)
					delta[count] = -25.0*k1[count]/108.0+125.0*k4[count]/108.0-65.0*k5[count]/27.0+125.0*k6[count]/54.0;

				rhs(m1, m2, 
				x0[id] + delta[0], 
				x1[id] + delta[1],
				x2[id] + delta[2], 
				x3[id] + delta[3], 
				x4[id] + delta[4], 
				x5[id] + delta[5], 
				x6[id] + delta[6], 
				x7[id] + delta[7],
				x8[id] + delta[8], 
				x9[id] + delta[9], 
				x10[id] + delta[10],
				&dx[0], &dx[1], &dx[2], &dx[3], &dx[4], &dx[5], &dx[6], &dx[7], &dx[8], &dx[9], &dx[10]);
			
				for(count=0;count<11;count++)
					k7[count] = h*dx[count];
				
				//Coefficinet k8:
				for(count=0;count<11;count++)
					delta[count] = 31.0*k1[count]/300.0+61.0*k5[count]/225.0-2.0*k6[count]/9.0+13.0*k7[count]/900.0;

				rhs(m1, m2, 
				x0[id] + delta[0], 
				x1[id] + delta[1],
				x2[id] + delta[2], 
				x3[id] + delta[3], 
				x4[id] + delta[4], 
				x5[id] + delta[5], 
				x6[id] + delta[6], 
				x7[id] + delta[7],
				x8[id] + delta[8], 
				x9[id] + delta[9], 
				x10[id] + delta[10],
				&dx[0], &dx[1], &dx[2], &dx[3], &dx[4], &dx[5], &dx[6], &dx[7], &dx[8], &dx[9], &dx[10]);
			
				for(count=0;count<11;count++)
					k8[count] = h*dx[count];


				//Coefficinet k9:
				for(count=0;count<11;count++)
					delta[count] = 2.0*k1[count]-53.0*k4[count]/6.0+704.0*k5[count]/45.0-107.0*k6[count]/9.0+67.0*k7[count]/90.0+3.0*k8[count];

				rhs(m1, m2, 
				x0[id] + delta[0], 
				x1[id] + delta[1],
				x2[id] + delta[2], 
				x3[id] + delta[3], 
				x4[id] + delta[4], 
				x5[id] + delta[5], 
				x6[id] + delta[6], 
				x7[id] + delta[7],
				x8[id] + delta[8], 
				x9[id] + delta[9], 
				x10[id] + delta[10],
				&dx[0], &dx[1], &dx[2], &dx[3], &dx[4], &dx[5], &dx[6], &dx[7], &dx[8], &dx[9], &dx[10]);
			
				for(count=0;count<11;count++)
					k9[count] = h*dx[count];

				//Coefficinet k10:
				for(count=0;count<11;count++)
					delta[count] = -91.0*k1[count]/108.0+23.0*k4[count]/108.0-976.0*k5[count]/135.0+311.0*k6[count]/54.0-19.0*k7[count]/60.0+17.0*k8[count]/6.0-k9[count]/12.0;

				rhs(m1, m2, 
				x0[id] + delta[0], 
				x1[id] + delta[1],
				x2[id] + delta[2], 
				x3[id] + delta[3], 
				x4[id] + delta[4], 
				x5[id] + delta[5], 
				x6[id] + delta[6], 
				x7[id] + delta[7],
				x8[id] + delta[8], 
				x9[id] + delta[9], 
				x10[id] + delta[10],
				&dx[0], &dx[1], &dx[2], &dx[3], &dx[4], &dx[5], &dx[6], &dx[7], &dx[8], &dx[9], &dx[10]);
			
				for(count=0;count<11;count++)
					k10[count] = h*dx[count];

				//Coefficinet k11:
				for(count=0;count<11;count++)
					delta[count] = 2383.0*k1[count]/4100.0-341.0*k4[count]/164.0+4496.0*k5[count]/1025.0-301.0*k6[count]/82.0+2133.0*k7[count]/4100.0+45.0*k8[count]/82.0+45.0*k9[count]/164.0+18.0*k10[count]/41.0;

				rhs(m1, m2, 
				x0[id] + delta[0], 
				x1[id] + delta[1],
				x2[id] + delta[2], 
				x3[id] + delta[3], 
				x4[id] + delta[4], 
				x5[id] + delta[5], 
				x6[id] + delta[6], 
				x7[id] + delta[7],
				x8[id] + delta[8], 
				x9[id] + delta[9], 
				x10[id] + delta[10],
				&dx[0], &dx[1], &dx[2], &dx[3], &dx[4], &dx[5], &dx[6], &dx[7], &dx[8], &dx[9], &dx[10]);
			
				for(count=0;count<11;count++)
					k11[count] = h*dx[count];

				//Coefficinet k12:
				for(count=0;count<11;count++)
					delta[count] = 3.0*k1[count]/205.0-6.0*k6[count]/41.0-3.0*k7[count]/205.0-3.0*k8[count]/41.0+3.0*k9[count]/41.0+6.0*k10[count]/41.0;

				rhs(m1, m2, 
				x0[id] + delta[0], 
				x1[id] + delta[1],
				x2[id] + delta[2], 
				x3[id] + delta[3], 
				x4[id] + delta[4], 
				x5[id] + delta[5], 
				x6[id] + delta[6], 
				x7[id] + delta[7],
				x8[id] + delta[8], 
				x9[id] + delta[9], 
				x10[id] + delta[10],
				&dx[0], &dx[1], &dx[2], &dx[3], &dx[4], &dx[5], &dx[6], &dx[7], &dx[8], &dx[9], &dx[10]);
			
				for(count=0;count<11;count++)
					k12[count] = h*dx[count];

				//Coefficinet k13:
				for(count=0;count<11;count++)
					delta[count] = -1777.0*k1[count]/4100.0-341.0*k4[count]/164.0+4496.0*k5[count]/1025.0-289.0*k6[count]/82.0+2193.0*k7[count]/4100.0+51.0*k8[count]/82.0+33.0*k9[count]/164.0+12.0*k10[count]/41.0+k12[count];

				rhs(m1, m2, 
				x0[id] + delta[0], 
				x1[id] + delta[1],
				x2[id] + delta[2], 
				x3[id] + delta[3], 
				x4[id] + delta[4], 
				x5[id] + delta[5], 
				x6[id] + delta[6], 
				x7[id] + delta[7],
				x8[id] + delta[8], 
				x9[id] + delta[9], 
				x10[id] + delta[10],
				&dx[0], &dx[1], &dx[2], &dx[3], &dx[4], &dx[5], &dx[6], &dx[7], &dx[8], &dx[9], &dx[10]);
			
				for(count=0;count<11;count++)
					k12[count] = h*dx[count];

				//Error estimation
				error=0;
				for(count=0;count<11;count++)
					error=pow(41.0*k1[count]/840.0+41.0*k11[count]/840.0-41.0*k12[count]/840.0-41.0*k13[count]/840.0,2)+error;
				
				if(sqrt(error)> tol)
				{
					h=h/2;
					end=2*end;
					goto start;
				}
				for(count=0;count<11;count++)
					delta[count] = 41.0*k1[count]/840.0+34.0*k6[count]/105.0+9.0*k7[count]/35.0+9.0*k8[count]/35.0+9.0*k9[count]/280.0+9.0*k10[count]/280.0+41.0*k11[count]/840.0;

				tau[id]=MainLoop*h*end;
				x0[id]=x0[id]+delta[0];
				x1[id]=x1[id]+delta[1];
				x2[id]=x2[id]+delta[2];
				x3[id]=x3[id]+delta[3];
				x4[id]=x4[id]+delta[4];
				x5[id]=x5[id]+delta[5];
				x6[id]=x6[id]+delta[6];
				x7[id]=x7[id]+delta[7];
				x8[id]=x8[id]+delta[8];
				x9[id]=x9[id]+delta[9];
				x10[id]=x10[id]+delta[10];

				loop++;
			}
		}
		else
		{
			return;
		}
	}
}

