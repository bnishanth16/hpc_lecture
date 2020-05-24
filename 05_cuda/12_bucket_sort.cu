// Name: Nishanth Baskaran
// Student ID: 19M15017
// HPSC Assignment-L5
#include <cstdio>
#include <cstdlib>
#include <vector>

__global__ void init(int *bucket) {
  int i= blockIdx.x * blockDim.x + threadIdx.x;
  bucket[i]=0;
}

__global__ void add(int *key,int *bucket){
  int i= blockIdx.x * blockDim.x + threadIdx.x;
  atomicAdd(&bucket[key[i]],1);
}

__global__ void sort(int *key,int *bucket){
  int i= blockIdx.x * blockDim.x + threadIdx.x;

  for (int j=0,k=0; k<=i; j++){
    key[i]=j;
    __syncthreads();
    k+=bucket[j];
    __syncthreads();
  } 
}

int main() {
  int n = 50;
  int range = 5;

  int *key, *bucket;
  cudaMallocManaged(&key,n*sizeof(int));
  cudaMallocManaged(&bucket,range*sizeof(int));

  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");

  //since the range and n is small, using only 1 block for parallelisation
  init<<<1,range>>>(bucket);
  add<<<1,n>>>(key,bucket);
  sort<<<1,n>>>(key,bucket);
  cudaDeviceSynchronize();
  
  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");

  cudaFree(key);
  cudaFree(bucket);
}
