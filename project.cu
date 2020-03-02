#include <iostream>

__global__
void benes(int *x,int *table, int n)
{
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  int idx = threadIdx.x;
  int a , b;
  int level = blockIdx.x+1;

  if(level==1){
    a = x[idx*2];
    b = x[idx*2+1];
  }
  else {
    a = x[table[idx*2 + (n*blockIdx.x)]];
    b = x[table[idx*2 + (n*blockIdx.x)]];
  }  
  
  x[idx*2+(level*n)] = a;
  x[idx*2+(level*n)+1] = b;


}

int main(void){
    int *x;
    int inputSize = 16;
    int routerPerCol = inputSize /2;
    int col = (log2(inputSize)*2-1);

    int *table ;
    
    cudaMallocManaged(&x, inputSize*col*sizeof(int));
    cudaMallocManaged(&table, inputSize*col*sizeof(int));
    for (int i = 0; i < inputSize; i++) {
        if (i==0)
            x[i] = 1;
        else
            x[i] = x[i-1]+1;
      }
      for (int i=0; i<col * routerPerCol; i++)
          table[i] = 0;
     cudaDeviceSynchronize();
     benes<<< 7, 8>>>(x, table, inputSize);
    // cudaDeviceSynchronize();
    
    cudaFree(x);
}
 
 