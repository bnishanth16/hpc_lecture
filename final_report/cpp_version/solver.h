// HPSC Final Report
// Name: Nishanth Baskaran
// Student id: 19M15017
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

void initialize(double *u, double *un, double *v, double *vn, double *p, double *pn, double *b)
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

void build_up_b(double *b, double *u, double *v)
{
    for (int i = 0; i < nx-1; i++)
    {
        for (int j = 0; j < ny-1; j++)
        {
            b[j*nx+i] = (rho*(1.0/dt*
                        ((u[j*nx+i+1]-u[j*nx+i-1])/(2*dx) + (v[(j+1)*nx+i]-v[(j-1)*nx+i])/(2*dy)) -
                        ((u[j*nx+i+1]-u[j*nx+i-1])/(2*dx))*((u[j*nx+i+1]-u[j*nx+i-1])/(2*dx))-
                        2*((u[(j+1)*nx+i]-u[(j-1)*nx+i])/(2*dy)*(v[j*nx+i+1]-v[j*nx+i-1])/(2*dx))-
                        ((v[(j+1)*nx+i]-v[(j-1)*nx+i])/(2*dy))*((v[(j+1)*nx+i]-v[(j-1)*nx+i])/(2*dy))));
        }
    }
}

void copy_function(double *final, double *initial)
{
    for (int i = 0; i < nx; i++)
    {
        for (int j = 0; j < ny; j++)
        {
            final[j*nx+i] = initial[j*nx+i];
        }
    } 
}

void pressure_poisson(double *p, double *pn, double *b)
{
    for (int iter = 0; iter < nit; iter++)
    {
        copy_function(pn,p);

        for (int i = 1; i < nx-1; i++)
        {
            for (int j = 1; j < ny-1; j++)
            {
                p[j*nx+i] = (((pn[j*nx+i+1]+pn[j*nx+i-1])*dy*dy + (pn[(j+1)*nx+i]+pn[(j-1)*nx+i])*dx*dx)/
                            (2*(dx*dx+dy*dy))-
                            dx*dx*dy*dy*b[j*nx+i]*rho/(2*(dx*dx+dy*dy)));
            }
        }

        // Boundary Condition
        for (int i = 0; i < nx; i++)
        {
            p[i] = p[nx+i]; 
            p[(nx-1)*nx+i] = 0.0; 
        }

        for (int j = 1; j < ny-1; j++)
        {
            p[j*nx] = p[j*nx+1]; 
            p[(j+1)*nx-1] = p[(j+1)*nx-2]; 
        }
    }
}

void cavity_flow(double *u,double *un, double *v, double *vn, double *p, double *pn, double *b)
{
    copy_function(un,u);
    copy_function(vn,v);
    build_up_b(b,u,v);
    pressure_poisson(p,pn,b);

    for (int i = 1; i < nx-1; i++)
    {
        for (int j = 1; j < ny-1; j++)
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
    }

    // Boundary Condition
    for (int i = 0; i < nx; i++)
    {
        u[i] = 0.0;     u[(nx-1)*nx+i] = 1.0;
        v[i] = 0.0;     v[(nx-1)*nx+i] = 0.0;
    }
    
    for (int j = 1; j < ny-1; j++)
    {
        u[j*nx] = 0.0;      u[(j+1)*nx-1] = 0.0;
        v[j*nx] = 0.0;      v[(j+1)*nx-1] = 0.0;
    }
    
}

void save_result(double *u, double *v, double *p) {
    ofstream file("results.txt");

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
 
