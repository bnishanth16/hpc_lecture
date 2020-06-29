# HPSC Final Report
# Name: Nishanth Baskaran
# Student id: 19M15017
import numpy
from matplotlib import pyplot, cm
from mpl_toolkits.mplot3d import Axes3D

with open("results.txt") as readf:
    ln = readf.readline()
    ny, nx = ln.split()
    ny= int(ny)
    nx= int(nx)
    size=ny*nx
    u=numpy.zeros(size)
    v=numpy.zeros(size)
    p=numpy.zeros(size)

    ln=readf.readline()
    k=0
    for j in ln.split():
        u[k]=float(j)
        k+=1

    ln=readf.readline()
    k=0
    for j in ln.split():
        v[k]=float(j)
        k+=1

    ln=readf.readline()
    k=0
    for j in ln.split():
        p[k]=float(j)
        k+=1
    
    u=u.reshape((ny,nx))
    v=v.reshape((ny,nx))
    p=p.reshape((ny,nx))
    X=numpy.linspace(0,2,nx)
    Y=numpy.linspace(0,2,ny)
    x,y=numpy.meshgrid(X,Y)
    fig=pyplot.figure(figsize=(11,7), dpi=100)
    pyplot.contourf(x,y,p,alpha=0.5,cmap=cm.viridis)
    pyplot.colorbar()
    pyplot.contour(x,y,p,cmap=cm.viridis)
    pyplot.quiver(x[::2,::2], y[::2,::2], u[::2,::2], v[::2,::2])
    pyplot.xlabel('X')
    pyplot.ylabel('Y')
    pyplot.savefig("cavity_flow.png")