import math
import numpy as np
from sympy import *
import matplotlib.pyplot as plt

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

#m is not important because it will cancel for our purposes, assume its "1" but actually can be any value and it won't change a thing
def kepler_graph(method):
    X=[1,0,0,1]
    ode = lambda t,x : kepler_ode(x)
    nsteps=100
    dt=2*np.pi/nsteps 
    if (method=='euler'):
        solution =np.array(euler_method(ode,X,dt,nsteps))
    if (method=='rk4'):
        solution =np.array(rk4_method(ode,X,dt,nsteps))
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
    print("{method_n} dE/E error is: {error}".format(method_n=method,error=E_error ))
    print("x final is {x} ,y final is {y}".format(x=xf,y=yf))
    print("vx final is {vx} ,vy final is {vy}".format(vx=vxf,vy=vyf))


def main():
    E_error_euler=kepler_graph('euler')
    E_error_kepler=kepler_graph('rk4')
    plt.grid()
    plt.legend()

main()

