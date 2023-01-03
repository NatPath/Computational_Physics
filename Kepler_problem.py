import math
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
            if (estimated_trunc_error==0):
                print(f'dt is {dt} while estimated_trunc_error is 0')
                break
            dt_old=dt
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

def euler_method(ode,initial_conditions,dt,nsteps):
    nsteps=nsteps+1
    solution= []
    solution.append(initial_conditions)
    for i in range(1,nsteps):
        prev_X = solution[i-1]
        prev_t= (i-1)*dt
        new_X = prev_X+dt*ode(prev_t,prev_X)
        solution.append(new_X)
    return np.array(solution)

def kepler_ode(X):
    x_dot=X[2]
    y_dot=X[3]
    x_ddot=-X[0]/((X[0]**2+X[1]**2)**(3/2))
    y_ddot=-X[1]/((X[0]**2+X[1]**2)**(3/2))
    return np.array([x_dot,y_dot,x_ddot,y_ddot])

def timing_decorator(func):
    def wrapper(*args,**kwargs):
        start = time.perf_counter()
        original_return_val = func(*args, **kwargs)
        end = time.perf_counter()
        print("time elapsed in ", func.__name__, ": ", end - start, sep='')
        return original_return_val
    return wrapper

#m is not important because it will cancel for our purposes, assume its "1" but actually can be any value and it won't change a thing
@timing_decorator
def kepler_graph(method,initial_conditions,nsteps,G=1,M=1):
    X=initial_conditions
    ode = lambda t,x : kepler_ode(x)
    v=np.sqrt(initial_conditions[2]**2+initial_conditions[3]**2)
    r=np.sqrt(initial_conditions[0]**2+initial_conditions[1]**2)
    orbital_energy=v**2/2-G*M/r
    a=G*M/(2*orbital_energy) #semi-major axis
    T=2*np.pi*np.sqrt(np.abs(a)**3/(G*M))
    print(f'Number of steps: {nsteps}')
    print(f'Total orbit time is {T}')
    dt=T/nsteps
    print(f'timestep is {dt}')
    if (method=='euler'):
        solution =np.array(euler_method(ode,X,dt,nsteps))
    if (method=='rk4'):
        solution =np.array(rk4_method(ode,X,dt,int(nsteps)))
    if (method=='rk4a'):
        di_threshold=0.0003
        solution =np.array(rk4a_method(ode,X,dt,T,di_threshold))
    t_arr = np.array([i*dt for i in range(nsteps+1)])
    plt.plot(solution.transpose()[0],solution.transpose()[1],marker='.',ls='--',label='numerical solution ' +method)
    x0=solution[0][0]
    vx0=solution[0][2]
    y0=solution[0][1]
    vy0=solution[0][3]
    xf=solution[-1][0]
    vxf=solution[-1][2]
    yf=solution[-1][1]
    vyf=solution[-1][3]
    E_kinetic_intitial=(vx0**2+vy0**2)/2
    V_intital=-(x0**2+y0**2)**(-1/2)
    E_initial=E_kinetic_intitial+V_intital
    E_kinetic_final=(vxf**2+vyf**2)/2
    V_final=-(xf**2+yf**2)**(-1/2)
    E_final=E_kinetic_final+V_final
    E_error=(E_final-E_initial)/E_initial
    print("Results:")
    print(f'{method} dE/E error is: {E_error}')
    print(f'x final is {xf} ,y final is {yf}')
    print(f'vx final is {vxf} ,vy final is {vyf}')
    return E_error

def hw2q3():
    initial_conditions_circular=[1,0,0,1]
    num_of_steps=100
    E_error_euler=kepler_graph('euler',initial_conditions_circular,num_of_steps)
    E_error_kepler=kepler_graph('rk4',initial_conditions_circular,num_of_steps)
    plt.grid()
    plt.legend()


def hw3q1():
    x0=1
    f=0.1
    vy0=1*f
    initial_conditions_circular=[x0,0,0,vy0]
    num_of_steps=26000
    E_error=kepler_graph('rk4a',initial_conditions_circular,num_of_steps)
    plt.grid()
    plt.legend()

def main():
    print('a run has started')
    hw3q1()


main()

