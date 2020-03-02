// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include <cuda.h>
#include <cuda_runtime.h>


//constucting 8x8 benes network
//four rows, 5 columns, 20 routers total


__global__
void benes(int N, int rows, int columns,  int* network){

}



int main(int argc, char *argv[]){
	int N = 16; //8x8 benes network

	int LUTsize = N*(log2((double)N)*2 - 2);

	int* LUT;
	cudaMallocManaged(&LUT,LUTsize*sizeof(int));
	
	int n;
	int M = N;
	int M2 = N;
	for (int i = 0; i < LUTsize/2; i+=N){
		//printf("i: %d\n",i);
		for (int j = 0; j < N; j += M2){
			M2 = M; 
			for (int k=0; k < M; k+=2){
				//printf("mem: %d\n",i+j+k);
				printf("n: %d\n",n);
				LUT[i+j+k] = n%N;
				LUT[i+j+k+1] = n%N + M/2;
				n++;
			}
			n = n*2;
			M = N/2;
		}
		
			
	}
	
	
	
	
	//benes<<<numBlocks,blockSize>>>(N, 4, 5, network);
	//cudaDeviceSynchronize();
	
	for (int i = 0; i < LUTsize/2; i++){
		if (i%N == 0) printf("\n");
		printf("%d ", LUT[i]);
	}
	printf("\n");
}