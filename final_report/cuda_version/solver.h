// HPSC Final Report
// Name: Nishanth Baskaran
// Student id: 19M15017
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

void initialize(double *u, double *un, double *v, double *vn, double *p, double *pn, double *b, int nx, int ny)
{
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            u[j * nx + i] = 0.0;    un[j * nx + i] = 0.0;
            v[j * nx + i] = 0.0;    vn[j * nx + i] = 0.0;
            p[j * nx + i] = 0.0;    pn[j * nx + i] = 0.0;
            b[j * nx + i] = 0.0;
        }
    }
}

__global__ void build_up_b(double *b, double *u, double *v, double rho, double dt, double dx, double dy, int nx, int ny)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;
    if (i>0 && i<nx-1 && j>0 && j<ny-1)
    {
        b[j*nx+i] = (rho*(1.0/dt*
                        ((u[j*nx+i+1]-u[j*nx+i-1])/(2*dx) + (v[(j+1)*nx+i]-v[(j-1)*nx+i])/(2*dy)) -
                        ((u[j*nx+i+1]-u[j*nx+i-1])/(2*dx))*((u[j*nx+i+1]-u[j*nx+i-1])/(2*dx))-
                        2*((u[(j+1)*nx+i]-u[(j-1)*nx+i])/(2*dy)*(v[j*nx+i+1]-v[j*nx+i-1])/(2*dx))-
                        ((v[(j+1)*nx+i]-v[(j-1)*nx+i])/(2*dy))*((v[(j+1)*nx+i]-v[(j-1)*nx+i])/(2*dy))));
    }
    __syncthreads();
}

void copy_function(double *final, double *initial, int nx, int ny)
{
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            final[j*nx+i] = initial[j*nx+i];
        }
    } 
}

__global__ void pressure_poisson(double *p, double *pn, double *b, double rho, double dt, double dx, double dy, int nx, int ny)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;

    if (i>0 && i<nx-1 && j>0 && j<ny-1)
    {
        p[j*nx+i] = (((pn[j*nx+i+1]+pn[j*nx+i-1])*dy*dy + (pn[(j+1)*nx+i]+pn[(j-1)*nx+i])*dx*dx)/
                            (2*(dx*dx+dy*dy))-
                            dx*dx*dy*dy*b[j*nx+i]*rho/(2*(dx*dx+dy*dy)));
    }

    __syncthreads();
}

__global__ void boundary_pressure(double *p, int nx, int ny)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;

    if (j==0 && i<nx)
    {
        p[j*nx+i] = p[(j+1)*nx+i];
    }
    else if (j==ny-1 && i<nx)
    {
        p[j*nx+i] = 0.0;
    }
    else if(i==nx-1 && j>0 && j<ny-1)
    {
        p[j*nx+i] = p[j*nx+i-1];
    }
    else if (i==0 && j>0 && j<ny-1)
    {
        p[j*nx+i] = p[j*nx+i+1];
    }
    
    __syncthreads();
}
__global__ void velocity_solver(double *u,double *un, double *v, double *vn, double *p, double *pn, double *b, double rho, double nu, double dt, double dx, double dy, int nx, int ny)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;

    if (i>0 && i<nx-1 && j>0 && j<ny-1)
    {
        u[j*nx+i] = (un[j*nx+i]-un[j*nx+i]*dt/dx*(un[j*nx+i]-un[j*nx+i-1])-
                         vn[j*nx+i]*dt/dy*(un[j*nx+i]-un[(j-1)*nx+i])-
                         dt/(2*rho*dx)*(p[j*nx+i+1]-p[j*nx+i-1])+
                         nu*(dt/(dx*dx)*(un[j*nx+i+1]-2*un[j*nx+i]+un[j*nx+i-1])+
                             dt/(dy*dy)*(un[(j+1)*nx+i]-2*un[j*nx+i]+un[(j-1)*nx+i])));

        v[j*nx+i] = (vn[j*nx+i]-un[j*nx+i]*dt/dx*(vn[j*nx+i]-vn[j*nx+i-1])-
                         vn[j*nx+i]*dt/dy*(vn[j*nx+i]-vn[(j-1)*nx+i])-
                         dt/(2*rho*dy)*(p[(j+1)*nx+i]-p[(j-1)*nx+i])+
                         nu*(dt/(dx*dx)*(vn[j*nx+i+1]-2*vn[j*nx+i]+vn[j*nx+i-1])+
                             dt/(dy*dy)*(vn[(j+1)*nx+i]-2*vn[j*nx+i]+vn[(j-1)*nx+i])));
    }
    __syncthreads();
}

__global__ void boundary_velocity (double *u, double *v,int nx, int ny)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;
    int j = threadIdx.y + blockIdx.y * blockDim.y;

    if (j==0 && i<nx)
    {
        u[j*nx+i] = 0.0;    v[j*nx+i] = 0.0;
    }
    else if (j==ny-1 && i<nx)
    {
        u[j*nx+i] = 1.0;    v[j*nx+i] = 0.0;
    }
    else if(i==nx-1 && j>0 && j<ny-1)
    {
        u[j*nx+i] = 0.0;    v[j*nx+i] = 0.0;
    }
    else if (i==0 && j>0 && j<ny-1)
    {
        u[j*nx+i] = 0.0;    v[j*nx+i] = 0.0;
    }
    __syncthreads();
}

void save_result(double *u, double *v, double *p, int nx, int ny) {
    ofstream file("results_cuda.txt");

    if(file.is_open()){
        file<<ny<<" "<<nx<<"\n";
        for (int j = 0; j < ny; j++)
        {
            for (int i = 0; i < nx; i++)
            {
                file<<u[j*nx+i]<<" ";
            }
        }
        file<<"\n";

        for (int j = 0; j < ny; j++)
        {
            for (int i = 0; i < nx; i++)
            {
                file<<v[j*nx+i]<<" ";
            }
        }
        file<<"\n";

        for (int j = 0; j < ny; j++)
        {
            for (int i = 0; i < nx; i++)
            {
                file<<p[j*nx+i]<<" ";
            }
        }
        file<<"\n";
    }

    file.close();
}
 