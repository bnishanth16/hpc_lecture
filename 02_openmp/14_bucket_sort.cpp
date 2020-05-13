// Name: Nishanth Baskaran
// Student ID: 19M15017
// HPSC Assignment-L2 
#include <cstdio>
#include <cstdlib>
#include <vector>

int main() {
  int n = 50;
  int range = 5;
  std::vector<int> key(n);
  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");

  std::vector<int> bucket(range); 
  
  #pragma omp parallel for //new code 
  for (int i=0; i<range; i++) {
    bucket[i] = 0;
  }
 
 #pragma omp parallel //new code
 #pragma omp sections firstprivate(bucket) //new code 
 {
   #pragma omp section //new code
  for (int i=0; i<n; i++) {
    bucket[key[i]]++;
  }
 
  #pragma omp section //new code
  for (int i=0, j=0; i<range; i++) { 
    for (; bucket[i]>0; bucket[i]--) {
      key[j++] = i;
    }
  }
 }
  

  #pragma omp parallel //new code
  for (int i=0; i<n; i++) {
    #pragma omp single // new code
    printf("%d ",key[i]);
  }
  printf("\n");
}
