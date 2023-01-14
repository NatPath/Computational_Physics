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
from mpl_toolkits.mplot3d import Axes3D

# Define the electric field vector at each point in 3D space
def plot_ElectricField(E,min,max,dim):
    if dim==3:
        X, Y, Z = np.meshgrid(np.linspace(min, max, 20), np.linspace(min, max, 20), np.linspace(min, max, 20))
        Ex, Ey, Ez = np.sin(X), np.cos(Y), np.sin(Z)

        # Set up a figure and 3D axis
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Plot the electric field vectors
        ax.quiver(X, Y, Z, Ex, Ey, Ez, length=0.5)

        # Show the plot
        plt.show()
    if dim==2:
        X, Y, Z = np.meshgrid(np.linspace(min, max, 20), np.linspace(min, max, 20), np.linspace(min, max, 20))
        Ex, Ey, Ez = np.sin(X), np.cos(Y), np.sin(Z)

        # Set up a figure and 3D axis
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        # Plot the electric field vectors
        ax.quiver(X, Y, Z, Ex, Ey, Ez, length=0.5)

        # Show the plot
        plt.show()

#Chat GPT did the following function
def plot_electric_field(electric_potential, surface):
  # Extract the electric field on the surface plane from the electric potential
  ex, ey, ez = np.gradient(-electric_potential)
  electric_field = [ex[surface], ey[surface], ez[surface]]

  # Calculate the magnitudes of the electric field vectors
  field_magnitude = np.linalg.norm(np.array([electric_field[0], electric_field[1]]), axis=0)

  # Normalize the field magnitudes to the range [0, 1]
  field_magnitude_normalized = (field_magnitude - field_magnitude.min()) / (field_magnitude.max() - field_magnitude.min())

  # Generate a colormap for the field magnitudes
  field_colormap = plt.cm.viridis(field_magnitude_normalized)

  # Define the x and y coordinates of the vectors
  x, y, z= np.meshgrid(np.arange(electric_field.shape[0]),
                     np.arange(electric_field.shape[1]),
                     np.arange(electric_field.shape[2]))

  # Set up the plot
  fig = plt.figure()
  ax = fig.add_subplot(111, projection='3d')

  # Plot the electric field vectors
  ax.quiver(x, y,z, electric_field[0], electric_field[1], electric_field[2], cmap="copper", pivot='mid')

  # Show the plot
  plt.show()



def plot4d(arr,N,title,min=-3.3,max=3.3,step=1.1):
    dx=(max-min)/N
    fig = plt.figure(figsize=(10,10),dpi=340)
    ax = fig.add_subplot(projection="3d")
    space = np.array([*product(range(N), range(N), range(N))]) # all possible triplets of numbers from 0 to N-1
    p=ax.scatter(space[:,0]*dx+min, space[:,1]*dx+min, space[:,2]*dx+min, c=arr,alpha=0.1,cmap="jet")
    #p=ax.scatter(np.arange(min,max,dx), np.arange(min,max,dx),np.arange(min,max,dx), c=arr,alpha=0.1,cmap="jet")
    fig.colorbar(p)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title(title)
    '''
    ax.set_xticks(np.arange(min, max, step))
    ax.set_yticks(np.arange(min, max, step))
    ax.set_zticks(np.arange(min, max, step))
    '''
    plt.show()

def plot3d(arr,N,title,min=-3.3,max=3.3,step=1.1):
    dx=(max-min)/N
    fig = plt.figure(figsize=(10,10),dpi=340)
    ax = fig.add_subplot()
    space = np.array([*product(range(N), range(N))]) # all possible triplets of numbers from 0 to N-1
    p=ax.scatter(space[:,0]*dx+min, space[:,1]*dx+min, c=arr,cmap="jet")
    #p=ax.scatter(np.arange(min,max,dx),np.arange(min,max,dx), c=arr,cmap="jet")
    fig.colorbar(p)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(title)
    '''
    ax.set_yticks(np.arange(min, max, step))
    ax.set_xticks(np.arange(min, max, step))
    '''
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
    '''
    if res:
        print('im on a box')
    '''
    return res

def in_holes(x,y,z,hole_half_width):
    count = 0
    if abs(x)<= hole_half_width:
        count+=1
    if abs(y)<= hole_half_width:
        count+=1
    if abs(z)<= hole_half_width:
        count+=1
    return count>=2

def q_2_box_potential(x,y,z):
    return 1

def q_3_box_potential(x,y,z):
    return 1/np.sqrt(x**2+y**2+z**2)

def grounded_sphere_potential_with_a_box_inside(x,y,z,R,h,f):
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
        f_res=f(x,y,z)
        return f_res, on_box(x,y,z,h)
        #return f(x,y,z), True
    return 0, False

def grounded_sphere_potential_with_a_box_with_a_hole_inside(x,y,z,R,h,f):
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
    if on_box(x,y,z,h) and not(in_holes(x,y,z,hole_half_width=1)):
        #print(f'im inside boxwith my homies x={x} y={y} z={z}')
        return f(x,y,z), True
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


def calculate_potential_modified(boundary,last_potential,rho,Nx,Ny,Nz,a=5,first_flag=False):
    dx=a/Nx
    dy=a/Ny
    dz=a/Nz
    x_range= np.linspace(0,a,Nx+1)
    y_range= np.linspace(0,a,Ny+1)
    z_range= np.linspace(0,a,Nz+1)

    #new_potential= np.array([[[boundary(x_range[i],y_range[j],z_range[k])[0] if boundary(x_range[i],y_range[j],z_range[k])[1] else average_around_tenzor(last_potential,np.array([i,j,k]))+(dx**2)*(1/6)*rho[(i,j,k)] for i in range(1,Nx-1)] for j in range(1,Ny-1)] for k in range(1,Nz-1)])
    new_potential=np.zeros((2*Nx,2*Ny,2*Nz))
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                boundary_res=boundary(x_range[i],y_range[j],z_range[k])
                if boundary_res[1] or first_flag:
                    #print(boundary_res)
                    new_potential[(i+Nx,j+Ny,k+Nz)]=boundary_res[0]
                else:
                    new_potential[(i+Nx,j+Ny,k+Nz)]= average_around_tenzor(last_potential,np.array([i+Nx,j+Ny,k+Nz]))+(dx**2)*(1/6)*rho(x_range[i],y_range[j],z_range[k])
                new_pot=new_potential[(i+Nx,j+Ny,k+Nz)]
                new_potential[(Nx-i-1,j+Ny,k+Nz)]=new_pot
                new_potential[(Nx+i,Ny-j-1,k+Nz)]=new_pot
                new_potential[(Nx-i-1,Ny-j-1,k+Nz)]=new_pot
                new_potential[(Nx+i,Ny+j,Nz-k-1)]=new_pot
                new_potential[(Nx-i-1,Ny+j,Nz-k-1)]=new_pot
                new_potential[(Nx+i,Ny-j-1,Nz-k-1)]=new_pot
                new_potential[(Nx-i-1,Ny-j-1,Nz-k-1)]=new_pot

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

def laplacian(phi,dx):
    dimz=phi.shape
    res=[[[(phi[(min(i+1,dimz[0]-1),j,k)]+phi[(max(i-1,0),j,k)]+phi[(i,min(j+1,dimz[1]-1),k)]+phi[(i,max(j-1,0),k)]+phi[(i,j,min(k+1,dimz[2]-1))]+phi[(i,j,max(k-1,0))]-6*phi[(i,j,k)])/(dx**2) for k in range(dimz[2])] for j in range(dimz[1])] for i in range(dimz[0])]
    return np.array(res)

def gradient(phi,dx):
    dimz=phi.shape
    print(f'dimz is {dimz}')
    res_x=np.array([[[(phi[(min(i+1,dimz[0]-1),j,k)]-phi[(i,j,k)])/dx for i in range(dimz[0])] for j in range(dimz[1])] for k in range(dimz[2])])
    res_y=np.array([[[(phi[(i,min(j+1,dimz[1]-1),k)]-phi[(i,j,k)])/dx for i in range(dimz[0])] for j in range(dimz[1])] for k in range(dimz[2])])
    res_z=np.array([[[(phi[(i,j,min(k+1,dimz[2]-1))]-phi[(i,j,k)])/dx for i in range(dimz[0])] for j in range(dimz[1])] for k in range(dimz[2])])
    print(f'res_x is {res_x}')
    return res_x,res_y,res_z

def hw5(N,q='2',initial=0):
    Nx,Ny,Nz=N,N,N
    half_width_of_box=1.5
    L=3.3
    R=sqrt(10)
    dx=L/N
    print(f'dx is {dx}')
    if q=='2' or q=='4':
        f=lambda x,y,z:1
        rho = lambda x,y,z:0
    if q=='3a':
        f = lambda x,y,z:1/np.sqrt(x**2+y**2+z**2) if (x,y,z)!=(0,0,0) else 1/np.sqrt(0.75*dx**2)
        rho= lambda x,y,z: 0
    if q=='3b':
        f = lambda x,y,z:1/np.sqrt(x**2+y**2+z**2) if (x,y,z)!=(0,0,0) else 1/np.sqrt(3*dx**2)
        rho= lambda x,y,z: 1/dx**3 if (x==0 and y==0 and z==0) else 0
    if q!='4': 
        boundary_potential= lambda x,y,z:grounded_sphere_potential_with_a_box_inside(x,y,z,R,half_width_of_box,f)
    else:
        boundary_potential= lambda x,y,z:grounded_sphere_potential_with_a_box_with_a_hole_inside(x,y,z,R,half_width_of_box,f)

    if type(initial)==int:
        potential=calculate_potential_modified(boundary_potential,np.zeros((2*N,2*N,2*N)),rho,N,N,N,a=L,first_flag=False)
    else:
        potential=calculate_potential_modified(boundary_potential,initial,rho,N,N,N,a=L,first_flag=False)
    i=0
    fraction_change=1e-8
    print(f'fraction change criterea is {fraction_change}')
    while i<200:
        print(i)
        i=i+1
        new_potential=calculate_potential_modified(boundary_potential,potential,rho,N,N,N,a=L)
        mse=np.mean((new_potential-potential)**2)
        print(f'mse is {mse}')
        if mse<fraction_change:
            print(f'broke at i={i}')
            potential=new_potential
            break
        potential=new_potential
    np.save(f'hw5_q{q}_potential_N{N}.npy',new_potential)
    mEx,mEy,mEz=np.gradient(potential,dx)
    Ez=-mEz

    rho=-1*laplacian(potential,dx)
    plot4d(potential,2*N,f'q{q} 3D potential visualization')
    plot3d(potential[N],2*N,f'q{q} Potential of x=0 cut')
    '''
    Plotting phi in 1d
    '''
    plt.figure(dpi=340)
    x_range=np.arange(-L,L,dx)
    plt.plot(x_range,potential[N][N],label=f'numerical result')
    if q=='3a' or q=='3b':
        half_range=np.linspace(-L,L,102)
        point_charge_theory=abs(half_range**(-1))
        plt.plot(half_range,point_charge_theory,label=f'potential of point charge in theory')
    plt.grid()
    plt.legend()
    plt.xlabel('z')
    plt.ylabel(r'$\phi$(x=0,y=0,z)')
    plt.title(f'q{q} Potential of x=0 y=0 cut')
    plt.show()


    plot4d(rho,2*N,f'q{q} 3D rho visualization')
    plot3d(rho[N+int(10*N/22)+1],2*N,rf'N+10*N/22+1 $\rho$(x=1.5,y,z) distribution in q{q} visualization')
    plot3d(rho[N+int(10*N/22)],2*N,rf'N+10*N/22 $\rho$(x=1.5,y,z) distribution in q{q} visualization')
    plot3d(rho[N+int(10*N/22)-1],2*N,rf'N+10*N/22-1 $\rho$(x=1.5,y,z) distribution in q{q} visualization')
    plot3d(rho[N],2*N,rf'$\rho$(x=0,y,z) distribution in q{q} visualization')
    '''
    Plotting rho in 2d
    '''
    plt.plot(x_range,rho[N+int(10*N/22)][N+int(10*N/22)])
    plt.grid()
    plt.xlabel('z')
    plt.ylabel(r'$\rho$(x=1.5,y=1.5,z)')
    plt.title(f'q{q} Charge distribution in x=1.5 y=1.5 cut')
    plt.show()

    plt.plot(x_range,rho[N+int(10*N/22)][N])
    plt.grid()
    plt.xlabel('z')
    plt.ylabel(r'$\rho$(x=1.5,y=0,z)')
    plt.title(f'q{q} Charge distribution in x=1.5 y=0 cut')
    plt.show()

    '''
    Plotting electring field in x=0 y=0 plane (the projection on z)
    '''
    plt.plot(x_range,Ez[N][N])
    plt.grid()
    plt.xlabel('z')
    plt.ylabel('Ez(x=0,y=0,z)')
    plt.title('Electric field at (x=0,y=0,z) projection on z axis')
    plt.show()
    return 

def q_check(q,use_prev_result):
    possible_Ns=np.array([22+11*i for i in range(100)])
    N = possible_Ns[2]

    if use_prev_result==True:
        initial=np.load(f'hw5_q{q}_potential_N{N}.npy')
        hw5(N,q=f'{q}',initial=initial)
    else:
        hw5(N,q=f'{q}')

def main():
    #q_check('2',True)
    #q_check('3a',False)
    #q_check('3b',True)
    q_check('4',False)


main()
