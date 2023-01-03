import math
from matplotlib import projections
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
import time
import timeit

def progressbar(it, prefix="", size=60, out=sys.stdout): # Python3.6+
    count = len(it)
    def show(j):
        x = int(size*j/count)
        print(f"{prefix}[{u'█'*x}{('.'*(size-x))}] {j}/{count}", end='\r', file=out, flush=True)
    show(0)
    for i, item in enumerate(it):
        yield item
        show(i+1)
    print("\n", flush=True, file=out)
def solve_heat_eq(f_boundary,Nx,dt,nsteps,k=1):
    solution=[]
    dx=1/Nx
    solution.append([f_boundary(i*dx,0) for i in range(Nx+1)])
    print(f'length of solution {len(solution)}')
    for i in range(nsteps):
    #for i in progressbar(range(nsteps)):
        solution.insert(len(solution),[solution[i][j]+k*dt*(1/(dx**2))*(solution[i][j+1]+solution[i][j-1]-2*solution[i][j]) for j in range(1,Nx)])
        # adding boundries
        solution[-1].insert(0,f_boundary(0,dt*(i+1)))
        solution[-1].insert(len(solution[-1]),f_boundary(1,dt*(i+1)))
    return np.array(solution,dtype='float')

def single_source_boundary(x,t,T0=1):
    if x==0:
        return T0
    if t==0 or x==1:
        return 0
    print('single source error - arguments are not part of boundary')

def hw3q3_point_source():
    T0=1
    Nx=100
    nsteps=10000
    dt_min=(1/(Nx**2))/10
    dt=dt_min
    solution= solve_heat_eq(single_source_boundary,Nx,dt,nsteps) 
    for i in range(len(solution)):
        if len(solution[i]) !=len(solution[0]):
            print(f'something wrong with setup {i}')
            print(f'len of i is {len(solution[i])} while len of 0 is {len(solution[0])}')
    print(solution)
    plt.figure(dpi=340)
    xmin,xmax,ymin,ymax=(0,1,0,dt*nsteps)
    plt.imshow(solution, cmap=plt.cm.jet,extent=[xmin,xmax,ymax,ymin], vmin=0, vmax=1)
    plt.xlabel('distance')
    plt.ylabel('time')


    plt.colorbar() 

def periodic_source_boundary(x,t,w=1,T0=100000000):
    if x==0:
        return T0*sin(w*t)
    if t==0 or x==1:
        return 0
    print('single source error - arguments are not part of boundary')

def hw3q3_periodic_source(w):
    T0=1
    Nx=100
    nsteps=100000
    dt_min=(1/(Nx**2))/10
    dt=dt_min
    print(f'dt is {dt}')
    pde = lambda x,t : periodic_source_boundary(x,t,w)
    solution= solve_heat_eq(pde,Nx,dt,nsteps) 
    for i in range(len(solution)):
        if len(solution[i]) !=len(solution[0]):
            print(f'something wrong with setup {i}')
            print(f'len of i is {len(solution[i])} while len of 0 is {len(solution[0])}')
    print(solution)
    plt.figure(dpi=340)
    xmin,xmax,ymin,ymax=(0,Nx*10,0,nsteps)
    plt.imshow(solution, cmap=plt.cm.jet,extent=[xmin,xmax,ymax,ymin], vmin=0, vmax=1)
    plt.xlabel('distance')
    plt.ylabel('time')
    plt.colorbar() 

def main():
    #hw3q3_point_source()
    hw3q3_periodic_source(0.1)
    hw3q3_periodic_source(0.00000001)

main()

