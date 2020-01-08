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
#define BLOCKDIM 32

__global__
void apply_w(int * data,int * result, double * filter, int n){
	int my_x=blockIdx.x*blockDim.x+threadIdx.x;
	int my_y=blockIdx.y*blockDim.y+threadIdx.y;
	int my_id=my_x*n+my_y;
	//If thread is outside of compute id threshhold it doesnt need to do anything
	if(my_x>=n||my_y>=n){
		return;
	}

	double sum=0;
	for(int i=0;i<5;i++){
		for(int j=0;j<5;j++){
			sum+=filter[i*5+j]*data[n*((n+(my_x+i-2))%n)+((n+(my_y+j-2))%n)];
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
		dim3 dimGrid(n/BLOCKDIM+1,n/BLOCKDIM+1);
		apply_w<<<dimGrid,dimBlock>>(dev_G,dev_res,dev_w,n);
		dev_temp=dev_res;
		dev_res=dev_G;
		dev_G=dev_temp;
	}

	//Bring results back to CPU Host
	cudaMemcpy(G,dev_G,n*n*sizeof(int),cudaMemcpyDefault);



}

int main(int argc, char* argv[]){

	int buffer[267289];
	FILE *ptr;
	double w[25]={0.004,0.016,0.026,0.016,0.004,0.016,0.071,0.117,0.071,0.016,0.026,0.117,0,0.117,0.026,0.016,0.071,0.117,0.071,0.016,0.004,0.016,0.026,0.016,0.004};
	ptr = fopen("conf-init.bin","rb");  // r for read, b for binary

	fread(buffer,sizeof(buffer),1,ptr); // read 10 bytes to our buffer
	fclose(ptr);

	ising(buffer,w,11,517);

	int test[267289];

        ptr = fopen("conf-11.bin","rb");  // r for read, b for binary

        fread(test,sizeof(test),1,ptr); // read 10 bytes to our buffer
        fclose(ptr);
	int a=0;
	for(int i=0;i<267289;i++){
		if(test[i]!=buffer[i]){
			a++;
		}
	}
	printf("Errors=%d\n",a);
}
