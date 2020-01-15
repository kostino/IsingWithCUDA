//! Ising model evolution
/*!

  \param G      Spins on the square lattice             [n-by-n]
  \param w      Weight matrix                           [5-by-5]
  \param k      Number of iterations                    [scalar]
  \param n      Number of lattice points per dim        [scalar]

  NOTE: Both matrices G and w are stored in row-major format.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <cuda.h>
#include "ising.h"
#define BLOCKDIM 16
#define THREADDIM 6 //This is (int)sqrt(MAXMOMENTSPERTHREAD)

__global__
void apply_w(int * data,int * result, double * filter, int n){

	__shared__ int G_sh[BLOCKDIM * THREADDIM+4][BLOCKDIM * THREADDIM+4];
	__shared__ double w_sh[5][5];

	if (threadIdx.x < 5 && threadIdx.y < 5)
			w_sh[threadIdx.x][threadIdx.y] = filter[5*threadIdx.x + threadIdx.y];

	for (int x = 0; x <THREADDIM+1; x ++)
			{
					for (int y = 0; y <THREADDIM+1; y ++)
					{
							if((threadIdx.x*(THREADDIM+1)+x)<(BLOCKDIM * THREADDIM+4)&&(threadIdx.y*(THREADDIM+1)+y)<(BLOCKDIM * THREADDIM+4))
							{
								G_sh[threadIdx.x*(THREADDIM+1)+x][threadIdx.y*(THREADDIM+1)+y] =data[n*((n+BLOCKDIM*THREADDIM*blockIdx.x + threadIdx.x*(THREADDIM+1)+x-2)%n)+(n+BLOCKDIM*THREADDIM*blockIdx.y + threadIdx.y*(THREADDIM+1)+y-2)%n];
							}
					}
			}
	__syncthreads();


	for(int x=0;x<THREADDIM;x++){
		for(int y=0;y<THREADDIM;y++){
			int my_x=(blockIdx.x*blockDim.x+threadIdx.x)*THREADDIM+x;
			int my_y=(blockIdx.y*blockDim.y+threadIdx.y)*THREADDIM+y;
			int my_id=my_x*n+my_y;
			//If thread is outside of compute id threshhold it doesnt need to do anything
			if(my_x>=n||my_y>=n){
				break;
			}
			double sum=0;
			for(int i=0;i<5;i++){
				for(int j=0;j<5;j++){
					sum+=w_sh[i][j]*G_sh[threadIdx.x*THREADDIM+x+i][threadIdx.y*THREADDIM+y+j`];
				}
			}
			if((sum<1e-5)&&(sum>-(1e-5))){
				result[my_id]=data[my_id];
			}
			else if(sum<0){
				result[my_id]=-1;
			}
			else{
				result[my_id]=1;
			}
		}
	}
}


void ising( int *G, double *w, int k, int n){
	int * dev_temp;
	int * dev_G;
	int * dev_res;
	double * dev_w;
	if(cudaMalloc(&dev_G,n*n*sizeof(int))!=cudaSuccess||cudaMalloc(&dev_res,n*n*sizeof(int))!=cudaSuccess||cudaMalloc(&dev_w,25*sizeof(double))!=cudaSuccess){
		printf("Error: could not allocate memory on device!");
		return;
	}
	//copy data to GPU Device
	cudaMemcpy(dev_G,G,n*n*sizeof(int),cudaMemcpyDefault);
	cudaMemcpy(dev_w,w,25*sizeof(double),cudaMemcpyDefault);


	//execute kernel
	for(int rep=0;rep<k;rep++){
		dim3 dimBlock(BLOCKDIM,BLOCKDIM);
		dim3 dimGrid(n/(BLOCKDIM*THREADDIM)+1,n/(BLOCKDIM*THREADDIM)+1);
		apply_w<<<dimGrid,dimBlock>>>(dev_G,dev_res,dev_w,n);
		dev_temp=dev_res;
		dev_res=dev_G;
		dev_G=dev_temp;
	}

	//Bring results back to CPU Host
	cudaMemcpy(G,dev_G,n*n*sizeof(int),cudaMemcpyDefault);



}
