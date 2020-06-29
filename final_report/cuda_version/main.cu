// HPSC Final Report
// Name: Nishanth Baskaran
// Student id: 19M15017
#include <iostream>
#include <vector>
#include <chrono>
#include "solver.h"

using namespace std;

int main()
{
  const int nx = 41;  const int ny = 41;
  const int nt = 700; const int nit = 50; 
  
  const double L = 2.0;
  const double dx = L/(nx-1);    
  const double dy = L/(ny-1);
  
  const double rho = 1.0;
  const double nu = 0.1;
  const double dt = 0.001;
    
  int size = nx*ny*sizeof(double);
  double *u, *un, *v, *vn, *p, *pn, *b;

  cudaMallocManaged(&u,size);   cudaMallocManaged(&un,size);
  cudaMallocManaged(&v,size);   cudaMallocManaged(&vn,size);
  cudaMallocManaged(&p,size);   cudaMallocManaged(&pn,size);
  cudaMallocManaged(&b,size);

  initialize(u,un,v,vn,p,pn,b,nx,ny);

  dim3 threadsPerBlock(128,1);
  dim3 blockNumber ((nx+threadsPerBlock.x - 1)/threadsPerBlock.x, (ny+threadsPerBlock.y - 1)/threadsPerBlock.y);

  auto t_initial = chrono::steady_clock::now();

  for (int iter = 0; iter < nt; iter++)
  {
    build_up_b<<<blockNumber,threadsPerBlock>>>(b,u,v,rho,dt,dx,dy,nx,ny);
    cudaDeviceSynchronize();

    for (int p_iter = 0; p_iter < nit; p_iter++)
    {
      pressure_poisson<<<blockNumber,threadsPerBlock>>>(p,pn,b,rho,dt,dx,dy,nx,ny);
      boundary_pressure<<<blockNumber,threadsPerBlock>>>(p,nx,ny);
      copy_function(pn,p,nx,ny);
      cudaDeviceSynchronize();
    }
    
    velocity_solver<<<blockNumber,threadsPerBlock>>>(u,un,v,vn,p,pn,b,rho,nu,dt,dx,dy,nx,ny);
    boundary_velocity<<<blockNumber,threadsPerBlock>>>(u,v,nx,ny);
    copy_function(un,u,nx,ny);    copy_function(vn,v,nx,ny);
    cudaDeviceSynchronize();
  }
  auto t_final = chrono::steady_clock::now();
  double time = chrono::duration<double>(t_final-t_initial).count();
  cout << "time = " << time << endl;
  save_result(u,v,p,nx,ny);

  cudaFree(u); cudaFree(un); cudaFree(v); cudaFree(vn); cudaFree(p); cudaFree(pn); cudaFree(b);
  return 0;
}

