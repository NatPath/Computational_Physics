import math
from matplotlib import projections
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
import time
import timeit

#rk4 capable of solving system of equations
def rk4_method(ode,initial_conditions,dt,nsteps):
    nsteps=nsteps+1
    solution= []
    solution.append(initial_conditions)
    for i in range(1,nsteps):
        prev_t=(i-1)*dt
        prev_X=solution[i-1]
        F1=ode(prev_t,prev_X)
        F2=ode(prev_t+dt/2,prev_X+dt*F1/2)
        F3=ode(prev_t+dt/2,prev_X+dt*F2/2)
        F4=ode(prev_t+dt,prev_X+dt*F3)
        new_X=prev_X+(dt/6)*(F1+2*F2+2*F3+F4)
        solution.append(new_X)
    return np.array(solution)
def rk4a_method(ode,initial_conditions,dt,T,di_threshold):
    solution= []
    solution.append(initial_conditions)
    i,t=(0,0)
    while t<T:
        max_number_of_attempts=200
        for j in range(max_number_of_attempts):
            single_step_res=rk4_method(ode,solution[i],dt,1)
            double_step_res=rk4_method(ode,solution[i],dt/2,2)
            estimated_trunc_error= np.linalg.norm(double_step_res[-1]-single_step_res[-1])
            dt_old=dt
            if (estimated_trunc_error==0):
                print(f'dt is {dt} while estimated_trunc_error is 0')
                break
            dt_est=dt*np.abs(di_threshold/estimated_trunc_error)**0.2
            S2=2
            S1=1/2
            if dt_est/S1>dt*S2:
                dt=S2*dt
            elif S1*dt_est<dt/S2:
                dt=dt/S2
            else:
                dt=S1*dt_est
            if estimated_trunc_error<di_threshold:
                print(f'interation {i} t= {t} broke at j=={j}')
                break
            if j==199:
                print('The max_number_of_attempts has been exceeded')
        t+=dt_old
        if t<T*1.05:
            solution.append(single_step_res[-1])
        i+=1
    return np.array(solution)

def gravity_central_force(planet_location,particle_location):
    DX= particle_location[0]-planet_location[0] #x distance from planet
    DY= particle_location[1]-planet_location[1] #y distance from planet
    x_ddot=-DX/((DX**2+DY**2)**(3/2))
    y_ddot=-DY/((DX**2+DY**2)**(3/2))
    return [x_ddot,y_ddot,0]

def particle_near_multiple_planets_ode(t,X,planets_motions):
    x_dot=X[3]
    y_dot=X[4]
    z_dot=X[5]
    x_ddot,y_ddot,z_ddot=np.sum([gravity_central_force(planet_motion(t),X) for planet_motion in planets_motions],axis=0)
    return np.array([x_dot,y_dot,z_dot,x_ddot,y_ddot,z_ddot])


def particle_near_planets_graph(initial_conditions,T,G=1,M=1):
    X=initial_conditions
    planets_motions=[lambda t:(np.sin(t),np.cos(t)),lambda t: (-np.sin(t),-np.cos(t))]
    ode = lambda t,X :particle_near_multiple_planets_ode(t,X,planets_motions)
    nsteps=1000
    dt=T/nsteps
    di_threshold=0.000001
    print(f'threshold is {di_threshold}')
    solution =np.array(rk4a_method(ode,X,dt,2*T,di_threshold))
    #t_arr = np.array([i*dt for i in range(nsteps+1)])
    plt.figure(dpi=240)
    ax=plt.axes(projection='3d')
    ax.scatter3D(solution.transpose()[0],solution.transpose()[1],solution.transpose()[2],marker='.',ls='--',label='numerical solution ' )
    plt.grid()
    plt.legend()
    plt.figure(dpi=240)
    plt.plot(solution.transpose()[0],solution.transpose()[1],marker='.',ls='--',label='numerical solution ' )
    plt.grid()
    plt.legend()

def hw3q2():
    initial_conditions=[0,0,0,-1,0,0]
    T=2*np.pi
    particle_near_planets_graph(initial_conditions,T)
    
def main():
    hw3q2()
main()