import math
from matplotlib import projections
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
import time
import timeit
import itertools
from itertools import product
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

def plot4d(arr,N,title):
    fig = plt.figure(figsize=(10,10),dpi=340)
    ax = fig.add_subplot(projection="3d")
    space = np.array([*product(range(N), range(N), range(N))]) # all possible triplets of numbers from 0 to N-1
    p=ax.scatter(space[:,0], space[:,1], space[:,2], c=arr,alpha=0.1,cmap="jet")
    fig.colorbar(p)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title(title)
    plt.show()

def plot3d(arr,N,title):
    fig = plt.figure(figsize=(10,10),dpi=340)
    ax = fig.add_subplot()
    space = np.array([*product(range(N), range(N))]) # all possible triplets of numbers from 0 to N-1
    p=ax.scatter(space[:,0], space[:,1], c=arr,cmap="jet")
    fig.colorbar(p)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(title)
    plt.show()

def inside_sphere(x,y,z,R):
    r_squared=x**2+y**2+z**2
    return r_squared<R**2

def outside_sphere(x,y,z,R):
    return not(inside_sphere(x,y,z,R))

def inside_box(x,y,z,h):
    return (abs(x)<=h and abs(y)<=h and abs(z)<=h)

def on_box(x,y,z,h):
    res=((abs(x)==h or abs(y)==h or abs(z)==h) and inside_box(x,y,z,h))
    return res

def q_2_box_potential(x,y,z):
    return 1

def q_3_box_potential(x,y,z):
    return 1/np.sqrt(x**2+y**2+z**2)

def grounded_sphere_potential_inside_box(x,y,z,R,h,f):
    '''
    gets coordinates, radius of sphere ,half_width of box and a function on the box
    returns
    1) the potential at the point in the first argument
    2) if the point belongs to the boundary
    '''
    '''
    if on_box(x,y,z,h):
        print('foo')
    '''
    if outside_sphere(x,y,z,R):
        #print(f'im outside sphere with my homies x={x} y={y} z={z}')
        return 0, True
    #if on_box(x,y,z,h):
    if inside_box(x,y,z,h):
        #print(f'im inside boxwith my homies x={x} y={y} z={z}')
        return f(x,y,z), on_box(x,y,z,h)
        #return f(x,y,z), True
    return 0, False

def legal_indices(tenzor,indices_around,indices):
    N=[tenzor.shape[i] for i in range(len(tenzor.shape))]
    if tuple(indices_around)==tuple(indices):
        return False
    for j in range(len(indices_around)):
        if indices_around[j]<0 or indices_around[j]>=N[j]:
            return False
    return True


def average_around_tenzor(tenzor,indices):
    N=[tenzor.shape[i] for i in range(len(tenzor.shape))]
    sum=0
    count=0
    for indices_delta in itertools.product(range(-1,2),repeat=len(N)):
        indices_around=indices+np.array(indices_delta)
        if legal_indices(tenzor,indices_around,indices):
            sum+=tenzor[tuple(indices_around)]
            count+=1
    return sum/count


def calculate_potential_modified(boundary,last_potential,rho,Nx,Ny,Nz,a=5):
    dx=a/(Nx-2)
    dy=a/(Ny-2)
    dz=a/(Nz-2)
    x_range= np.linspace(-dx,a,Nx)
    y_range= np.linspace(-dy,a,Ny)
    z_range= np.linspace(-dz,a,Nz)

    #new_potential= np.array([[[boundary(x_range[i],y_range[j],z_range[k])[0] if boundary(x_range[i],y_range[j],z_range[k])[1] else average_around_tenzor(last_potential,np.array([i,j,k]))+(dx**2)*(1/6)*rho[(i,j,k)] for i in range(1,Nx-1)] for j in range(1,Ny-1)] for k in range(1,Nz-1)])
    new_potential=np.zeros((2*Nx,2*Ny,2*Nz))
    first_flag=True
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                boundary_res=boundary(x_range[i],y_range[j],z_range[k])
                if boundary_res[1] or first_flag:
                    #print(boundary_res)
                    new_potential[(i+Nx,j+Ny,k+Nz)]=boundary_res[0]
                    first_flag=False
                else:
                    new_potential[(i,j,k)]= average_around_tenzor(last_potential,np.array([i+Nx,j+Ny,k+Nz]))+(dx**2)*(1/6)*rho(x_range[i],y_range[j],z_range[k])
    #fix artificial boundary layer to fit the symetry 
    '''
    for i in range(2):
        for j in range(2):
            for k in range(2):
                new_potential[(i,j,k)]=new_potential[(2-i,2-j,2-k)]
    '''

    return new_potential
def calculate_potential(boundary,last_potential,rho,Nx,Ny,Nz,a=5):
    dx=a/(Nx-2)
    dy=a/(Ny-2)
    dz=a/(Nz-2)
    x_range= np.linspace(-dx,a,Nx)
    y_range= np.linspace(-dy,a,Ny)
    z_range= np.linspace(-dz,a,Nz)

    #new_potential= np.array([[[boundary(x_range[i],y_range[j],z_range[k])[0] if boundary(x_range[i],y_range[j],z_range[k])[1] else average_around_tenzor(last_potential,np.array([i,j,k]))+(dx**2)*(1/6)*rho[(i,j,k)] for i in range(1,Nx-1)] for j in range(1,Ny-1)] for k in range(1,Nz-1)])
    new_potential=np.zeros((Nx,Ny,Nz))
    first_flag=True
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                boundary_res=boundary(x_range[i],y_range[j],z_range[k])
                if boundary_res[1] or first_flag:
                    #print(boundary_res)
                    new_potential[(i,j,k)]=boundary_res[0]
                    first_flag=False
                else:
                    new_potential[(i,j,k)]= average_around_tenzor(last_potential,np.array([i,j,k]))+(dx**2)*(1/6)*rho(x_range[i],y_range[j],z_range[k])
    #fix artificial boundary layer to fit the symetry 
    '''
    for i in range(2):
        for j in range(2):
            for k in range(2):
                new_potential[(i,j,k)]=new_potential[(2-i,2-j,2-k)]
    '''

    return new_potential


def hw5(N,q='2'):
    Nx,Ny,Nz=N,N,N
    half_width_of_box=1.5
    L=3.3
    R=sqrt(10)
    dx=L/(N-2)
    print(f'dx is {dx}')
    if q=='2':
        f=lambda x,y,z:1
        rho = lambda x,y,z:0
    if q=='3a':
        f = lambda x,y,z:1/np.sqrt(x**2+y**2+z**2)
        rho= lambda x,y,z: 0
    if q=='3b':
        f = lambda x,y,z:1/np.sqrt(x**2+y**2+z**2)
        rho= lambda x,y,z: 1 if (x==0 and y==0 and z==0) else 0
    
    boundary_potential= lambda x,y,z:grounded_sphere_potential_inside_box(x,y,z,R,half_width_of_box,f)
    potential=calculate_potential(boundary_potential,np.zeros((N,N,N)),rho,N,N,N,L)
    i=0
    fraction_change=0.00001
    while i<200:
        print(i)
        i=i+1
        new_potential=calculate_potential(boundary_potential,potential,rho,N,N,N,L)
        if np.mean((new_potential-potential)**2)<fraction_change:
            print(f'broke at i={i}')
            break
        potential=new_potential
    np.save(f'hw5_q{q}_potential.npy',new_potential)
    E=-1*np.gradient(potential)
    plot4d(potential,N,f'q{q} 3D potential visualization')
    plot3d(potential[1],N,f'q{q} Potential of x=0 cut')
    return 


def main():
    possible_Ns=np.array([24+11*i for i in range(100)])
    N = possible_Ns[0]
    hw5(N,q='2')
    arr=np.load('hw5_q2_potential.npy')

    return arr


arr= main()
