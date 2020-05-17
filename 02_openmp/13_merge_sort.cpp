// Name: Nishanth Baskaran
// Student ID: 19M15017
// HPSC Assignment-L2 
#include <cstdio>
#include <cstdlib>
#include <vector>

template<class T>
void merge(std::vector<T>& vec, int begin, int mid, int end) {
  std::vector<T> tmp(end-begin+1);
  int left = begin;
  int right = mid+1;
  
    for (int i=0; i<tmp.size(); i++) {
    if (left > mid)
      tmp[i] = vec[right++];
    else if (right > end)
      tmp[i] = vec[left++];
    else if (vec[left] <= vec[right])
      tmp[i] = vec[left++];
    else
      tmp[i] = vec[right++];     
  }
  
  for (int i=0; i<tmp.size(); i++) {
    vec[begin++] = tmp[i];
  }
}

template<class T>
void merge_sort(std::vector<T>& vec, int begin, int end) {
  if(begin < end) {
    int mid = (begin + end) / 2;

    #pragma omp task shared(vec) //New code
    merge_sort(vec, begin, mid);

    #pragma omp task shared(vec) //New code
    merge_sort(vec, mid+1, end);

    #pragma omp taskwait //New code
    merge(vec, begin, mid, end);
  }
}

int main() {
  int n = 20;
  std::vector<int> vec(n);
  
  for (int i=0; i<n; i++) {
    vec[i] = rand() % (10 * n);
    printf("%d ",vec[i]);
  }
  printf("\n");

  #pragma omp parallel // new code
  {
  #pragma omp single // new code
  merge_sort(vec, 0, n-1);
  }

  #pragma omp parallel //New code
  for (int i=0; i<n; i++) {
    #pragma omp single //New code
    printf("%d ",vec[i]);
  }
  printf("\n");
}
