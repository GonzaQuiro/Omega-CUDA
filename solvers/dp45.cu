// Version 1.0.0 CUDA-C: Runge-Kutta 4 method
// Dr. Gonzalo Damián Quiroga
#include "../ODEs/rhs.h"

void dp45(double m1, double m2, double *tau, double *x0, double *x1, double *x2, double *x3, double *x4, double *x5, double *x6, double *x7, double *x8, double *x9, double *x10, int n_ic, double h, double tol, int MainLoop)
{
	const int id = blockIdx.x*blockDim.x + threadIdx.x;
	double k1[11], k2[11], k3[11], k4[11], k5[11], k6[11], k7[11];
	double dx[11], delta[11];
	double error;
	int count;
	int end=1;
	int loop=0;

	if (id < n_ic)
	{	
		//Stop Condition
		if(x0[id]<0.05)
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
					delta[count] = k1[count]/5.0;

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
					delta[count] = 3.0*k1[count]/40.0 + 9.0*k2[count]/40.0;

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
					delta[count] = 44.0*k1[count]/45.0-56.0*k2[count]/15.0+32.0*k3[count]/9.0;

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
					delta[count] = 19372.0*k1[count]/6561.0-25360.0*k2[count]/2187.0+64448.0*k3[count]/6561.0-212.0*k4[count]/729.0;

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
					delta[count] = 9017.0*k1[count]/3168.0-355.0*k2[count]/33.0+46732.0*k3[count]/5247.0+49.0*k4[count]/176.0-5103.0*k5[count]/18656.0;

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
					delta[count] = 35.0*k1[count]/384.0+500.0*k3[count]/1113.0+125.0*k4[count]/192.0-2187.0*k5[count]/6784.0+11.0*k6[count]/84.0;

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
				
				//Error estimation
				error=0;
				for(count=0;count<11;count++)
					error=pow(-71.0*k1[count]/57600.0+71.0*k3[count]/16695.0-71.0*k4[count]/1920.0+17253.0*k5[count]/339200.0-22.0*k6[count]/525.0+k7[count]/40.0,2)+error;
				
				if(sqrt(error)> tol)
				{
					h=h/2;
					end=2*end;
					goto start;
				}
				
				for(count=0;count<11;count++)
					delta[count] = 35.0*k1[count]/384.0+500.0*k3[count]/1113.0+125.0*k4[count]/192.0-2187.0*k5[count]/6784.0+11.0*k6[count]/84.0;

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

