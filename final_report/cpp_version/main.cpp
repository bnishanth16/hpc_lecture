// HPSC Final Report
// Name: Nishanth Baskaran
// Student id: 19M15017
#include "parameter.h"
#include "solver.h"
#include <vector>
#include <chrono>

using namespace std;
int main()
{
  int size = nx*ny*sizeof(double);
  double *u= (double*)malloc(size);       double *un= (double*)malloc(size);
  double *v= (double*)malloc(size);       double *vn= (double*)malloc(size);
  double *p= (double*)malloc(size);       double *pn= (double*)malloc(size);
  double *b= (double*)malloc(size);

  initialize(u,un,v,vn,p,pn,b);

  auto t_initial = chrono::steady_clock::now();
  for (int i = 0; i < nt; i++)
  {
      cavity_flow(u,un,v,vn,p,pn,b);
  }
  auto t_final = chrono::steady_clock::now();
  double time = chrono::duration<double>(t_final-t_initial).count();
  cout << "time = " << time << endl;
  save_result(u,v,p);
  free(u); free(un); free(v); free(vn); free(p); free(pn); free(b);
  return 0;
}

