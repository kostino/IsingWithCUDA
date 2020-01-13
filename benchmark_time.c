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
#include <time.h>

#include "ising.h"


int main(int argc, char* argv[]){
    
    double start,end;
    
    int n=2;
    
	int buffer[267289];
    
    int data[1069156];
    
    
	FILE *ptr;
	double w[25]={0.004,0.016,0.026,0.016,0.004,0.016,0.071,0.117,0.071,0.016,0.026,0.117,0,0.117,0.026,0.016,0.071,0.117,0.071,0.016,0.004,0.016,0.026,0.016,0.004};
	ptr = fopen("conf-init.bin","rb");  // r for read, b for binary

	fread(buffer,sizeof(buffer),1,ptr); // read 10 bytes to our buffer
	fclose(ptr);
    
    for(int i=0;i<n*n*517*517;i++){
        data[i]=buffer[i%(517*517)];
    }
    
    start=clock();
    
	ising(data,w,1,517*n);

    end=clock();
    printf("time= %f sec\n",((double) (end - start)) / CLOCKS_PER_SEC);
}
