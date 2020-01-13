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
void ising( int *G, double *w, int k, int n){
  double sum=0;
  int *temp;
  int changes=0;
  int *res=(int*)calloc(n*n,sizeof(double));
  for(int reps=0;reps<k;reps++){
    changes=0;
  for(int x=0;x<n;x++){
    for(int y=0;y<n;y++){
        sum=0;
        for(int i=0;i<5;i++){
            for(int j=0;j<5;j++){
                sum+=w[i*5+j]*G[((n+(x+i-2))%n)*n+(n+(y+j-2))%n];
            }
        }
        if(sum==0||fabs(sum)<1e-5){
            res[x*n+y]=G[x*n+y];
            
        }
        else if(sum<0){
            res[x*n+y]=-1;
            
        }
        else{
            res[x*n+y]=1;
            
        }
        changes+=(res[x*n+y]!=G[x*n+y]);
    }
  }
      temp=G;
      G=res;
      res=temp;
  
  if(changes==0){
      break;
  }
  }
  if(k%2==1){
  memcpy(G,res,n*n*sizeof(int));
  }
 
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
