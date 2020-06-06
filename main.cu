// Version 1.0.0 CUDA-C: Omega
// Dr. Gonzalo Damián Quiroga
// Universidad Nacional de Córdoba

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <stdio.h>
#include "common/handle.h"
#include "common/menu.h"
#include "kernel/kernel.h"
#include "kernel/kernel.cu"
#include "common/menu.c"

struct record
{
   double t,x0,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11;
};

int main(int argc,char **argv)
{	
	//Start Program
	menu();
	//Setting the Initial conditions
	int n_ic;
	printf("Set the number of initial conditions: ");
	scanf("%d", &n_ic);

	//Define parameters
	double m1,m2;	
	printf("Set the black holes masses, m1 m2: ");
	scanf("%lf %lf", &m1, &m2);

	//Solver settings
	double h;
	double tol=1e-10; //tolerance for solver
	double final_time; // should be a high value;
	
	//Integration steps
	printf("Final Time: ");
	scanf("%lf", &final_time);
	printf("Initial step size: ");
	scanf("%lf", &h);

	int MFlag;
	//Flag Solver
	printf("Select the numerical method [1]-RK4 [2]-DP45 [3]-RKF78: ");
	scanf("%d", &MFlag);

	// Host input/output vectors	
	double *h_x0, *h_x1, *h_x2, *h_x3,*h_x4, *h_x5,*h_x6, *h_x7, *h_x8,*h_x9, *h_x10, *h_tau;
	 
	// Device input/output vectors  
	double *d_x0, *d_x1, *d_x2, *d_x3,*d_x4, *d_x5,*d_x6, *d_x7, *d_x8,*d_x9, *d_x10, *d_tau;
	
	// Size, in bytes, of each vector
	double nBytes = n_ic*sizeof(double);

	// Allocate memory for each vector on host
	h_tau= (double *)malloc(nBytes);
	h_x0= (double *)malloc(nBytes);
	h_x1= (double *)malloc(nBytes);
	h_x2= (double *)malloc(nBytes);
	h_x3= (double *)malloc(nBytes);
	h_x4= (double *)malloc(nBytes);
	h_x5= (double *)malloc(nBytes);
	h_x6= (double *)malloc(nBytes);
	h_x7= (double *)malloc(nBytes);
	h_x8= (double *)malloc(nBytes);
	h_x9= (double *)malloc(nBytes);
	h_x10= (double *)malloc(nBytes);
	
	// Allocate memory for each vector on GPU
	printf("Allocating device memory on host..\n");
	HANDLE_ERROR(cudaMalloc((void **)&d_tau,nBytes));
	HANDLE_ERROR(cudaMalloc((void **)&d_x0,nBytes));
	HANDLE_ERROR(cudaMalloc((void **)&d_x1,nBytes));
	HANDLE_ERROR(cudaMalloc((void **)&d_x2,nBytes));
	HANDLE_ERROR(cudaMalloc((void **)&d_x3,nBytes));
	HANDLE_ERROR(cudaMalloc((void **)&d_x4,nBytes));
	HANDLE_ERROR(cudaMalloc((void **)&d_x5,nBytes));
	HANDLE_ERROR(cudaMalloc((void **)&d_x6,nBytes));
	HANDLE_ERROR(cudaMalloc((void **)&d_x7,nBytes));
	HANDLE_ERROR(cudaMalloc((void **)&d_x8,nBytes));
	HANDLE_ERROR(cudaMalloc((void **)&d_x9,nBytes));
	HANDLE_ERROR(cudaMalloc((void **)&d_x10,nBytes));

	//Allocate memory of struct data
	struct record data;
	
	// Initial conditions on host
	FILE *fic = fopen("ic/ic.txt", "r");
	if (fic == NULL)
	{
		perror("Error: can't open ic.txt.");
		return -1;
	}

	//Read the initial conditions
	int count;
	for (count = 0; count < n_ic; count++)
	{
		fscanf(fic, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &h_x0[count], &h_x1[count], &h_x2[count], &h_x3[count], &h_x4[count], &h_x5[count], &h_x6[count], &h_x7[count], &h_x8[count], &h_x9[count], &h_x10[count]);
	}
	//Close the file
	fclose(fic);
	//Set initial time
	for (count = 0; count < n_ic; count++)
	{
		h_tau[count]=0.0;
	}

	//Set the Block and grid Size
	int blockSize, gridSize,minGridSize; 

	// Number of threads in each thread block.
	//Suggested block size to achieve maximum occupancy.
	if(MFlag==1)
		cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, rk4, 0, n_ic);
	else if(MFlag==2)
		cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, dp45, 0, n_ic);
	else if(MFlag==3)
		cudaOccupancyMaxPotentialBlockSize(&minGridSize, &blockSize, rkf78, 0, n_ic);

	gridSize = (n_ic + blockSize - 1) / blockSize;
	dim3 dimBlock(blockSize,1,1);
	dim3 dimGrid(gridSize,1,1);
	
	// Copy host vectors to device 
	printf("Copying to device..\n");
	HANDLE_ERROR(cudaMemcpy(d_tau,h_tau,nBytes,cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_x0,h_x0,nBytes,cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_x1,h_x1,nBytes,cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_x2,h_x2,nBytes,cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_x3,h_x3,nBytes,cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_x4,h_x4,nBytes,cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_x5,h_x5,nBytes,cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_x6,h_x6,nBytes,cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_x7,h_x7,nBytes,cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_x8,h_x8,nBytes,cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_x9,h_x9,nBytes,cudaMemcpyHostToDevice));
	HANDLE_ERROR(cudaMemcpy(d_x10,h_x10,nBytes,cudaMemcpyHostToDevice));

	//Evol: Time Start
   	clock_t start_d=clock();
	
	//Open fout
	//Configuration
	int format;
	printf("Set output format [1]-text [2]-binary: ");
	scanf("%d", &format);

	FILE *fout;
	if(format==1)
		fout = fopen("output/output.txt", "w");
	else if(format==2)
		fout = fopen("output/output.bin", "wb");


	if (fout == NULL)
	{
		perror("Error: the output folder does not exist.");
		return -1;
	}
	else
	{
		printf("Output file created.. \n");
	}

	
	printf("Evolving systems.. \n");

	int loop, loopMax=(final_time/h)+1;
	for(loop=1;loop<=loopMax;loop++)
	{
		printf("\r%.3lf / %.3lf", loop*h, final_time);
		//Print the initial conditions in the output file
		if(format==1)
		{	
			for (count = 0; count < n_ic; count++)
			{
				fprintf(fout, "%.8lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n", h_tau[count], h_x0[count], h_x1[count], h_x2[count], h_x3[count], h_x4[count], h_x5[count], h_x6[count], h_x7[count], h_x8[count], h_x9[count], h_x10[count]);	
			}
		}
		else if(format==2)
		{	
			for (count = 0; count < n_ic; count++)
			{
				data.t =h_tau[count];
				data.x0= h_x0[count]; 
				data.x1=h_x1[count];
				data.x2= h_x2[count]; 
				data.x3=h_x3[count];
				data.x4= h_x4[count];
				data.x5= h_x5[count];
				data.x6= h_x6[count];
				data.x7= h_x7[count];
				data.x8= h_x8[count];
				data.x9= h_x9[count];
				data.x10= h_x10[count];	
				fwrite(&data, sizeof(struct record), 1, fout); 
			}
		}

		// Executing kernel
		if(MFlag==1)
			rk4<<<gridSize,blockSize>>>(m1, m2, d_tau, d_x0,d_x1,d_x2,d_x3,d_x4,d_x5,d_x6,d_x7,d_x8,d_x9,d_x10, n_ic, h, loop);
		else if(MFlag==2)
			dp45<<<gridSize,blockSize>>>(m1, m2, d_tau, d_x0,d_x1,d_x2,d_x3,d_x4,d_x5,d_x6,d_x7,d_x8,d_x9,d_x10, n_ic, h, tol, loop);	
		else if(MFlag==3)
			rkf78<<<gridSize,blockSize>>>(m1, m2, d_tau, d_x0,d_x1,d_x2,d_x3,d_x4,d_x5,d_x6,d_x7,d_x8,d_x9,d_x10, n_ic, h, tol, loop);	
		//cudaThreadSynchronize();
		cudaDeviceSynchronize();

		// Copy array back to host	
		HANDLE_ERROR(cudaMemcpy(h_tau,d_tau,nBytes,cudaMemcpyDeviceToHost));
		HANDLE_ERROR(cudaMemcpy(h_x0,d_x0,nBytes,cudaMemcpyDeviceToHost));
		HANDLE_ERROR(cudaMemcpy(h_x1,d_x1,nBytes,cudaMemcpyDeviceToHost));
		HANDLE_ERROR(cudaMemcpy(h_x2,d_x2,nBytes,cudaMemcpyDeviceToHost));
		HANDLE_ERROR(cudaMemcpy(h_x3,d_x3,nBytes,cudaMemcpyDeviceToHost));
		HANDLE_ERROR(cudaMemcpy(h_x4,d_x4,nBytes,cudaMemcpyDeviceToHost));
		HANDLE_ERROR(cudaMemcpy(h_x5,d_x5,nBytes,cudaMemcpyDeviceToHost));
		HANDLE_ERROR(cudaMemcpy(h_x6,d_x6,nBytes,cudaMemcpyDeviceToHost));
		HANDLE_ERROR(cudaMemcpy(h_x7,d_x7,nBytes,cudaMemcpyDeviceToHost));
		HANDLE_ERROR(cudaMemcpy(h_x8,d_x8,nBytes,cudaMemcpyDeviceToHost));
		HANDLE_ERROR(cudaMemcpy(h_x9,d_x9,nBytes,cudaMemcpyDeviceToHost));
		HANDLE_ERROR(cudaMemcpy(h_x10,d_x10,nBytes,cudaMemcpyDeviceToHost));
		
	}
	printf("\rEvol Done!                                                  \n");
	
	//Evol: Time Ends
   	clock_t end_d = clock();
   	double time_spent = (double)(end_d-start_d)/CLOCKS_PER_SEC;

	//Close the output file
	fclose(fout);

	// Release device memory
   	cudaFree(d_x0);
	cudaFree(d_x1);
	cudaFree(d_x2);
	cudaFree(d_x3);
	cudaFree(d_x4);
	cudaFree(d_x5);
	cudaFree(d_x6);
	cudaFree(d_x7);
	cudaFree(d_x8);
	cudaFree(d_x9);
	cudaFree(d_x10);

	// Release host memory
	free(h_x0);	
	free(h_x1);
	free(h_x2);
	free(h_x3);
	free(h_x4);	
	free(h_x5);
	free(h_x6);
	free(h_x7);
	free(h_x8);	
	free(h_x9);
	free(h_x10);

	//Wall Time
	printf("Runtime: %f sec \n", time_spent);
		
	return 0;
}
