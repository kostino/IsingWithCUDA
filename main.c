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
    
	int buffer[267289];
	FILE *ptr;
	double w[25]={0.004,0.016,0.026,0.016,0.004,0.016,0.071,0.117,0.071,0.016,0.026,0.117,0,0.117,0.026,0.016,0.071,0.117,0.071,0.016,0.004,0.016,0.026,0.016,0.004};
	ptr = fopen("conf-init.bin","rb");  // r for read, b for binary

	fread(buffer,sizeof(buffer),1,ptr); // read 10 bytes to our buffer
	fclose(ptr);
    
    start=clock();
    
	ising(buffer,w,11,517);

    end=clock();
    
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
    printf("time= %f sec\n",((double) (end - start)) / CLOCKS_PER_SEC);
}
