import math
import numpy as np
from sympy import *
import matplotlib.pyplot as plt

def taylor(delta,vals,N):
    monoms = [ (delta**n)*vals[n]/math.factorial(n) for n in range(N)]
    return sum(monoms)

def euler_method(ode,initial_conditions,dt,nsteps):
    nsteps=nsteps+1
    n= len(initial_conditions)
    solution= [[] for i in range(nsteps)]
    solution[0]=[*initial_conditions,ode(0,*initial_conditions)]
    for i in range(1,nsteps):
        for j in range(1,n+1):
            solution[i].insert(0,taylor(dt,solution[i-1][n-j:],2))
        solution[i].append(ode(i*dt,*solution[i]))
    return np.array(solution)

def super_euler_method(ode,initial_conditions,dt,nsteps):
    nsteps=nsteps+1
    n= len(initial_conditions)
    solution= [[] for i in range(nsteps)]
    solution[0]=[*initial_conditions,ode(0,*initial_conditions)]
    for i in range(1,nsteps):
        for j in range(1,n+1):
            solution[i].insert(0,taylor(dt,solution[i-1][n-j:],j+1))
        solution[i].append(ode(i*dt,*solution[i]))
    return np.array(solution)

def rk2_method(ode,initial_conditions,dt,nsteps):
    nsteps=nsteps+1
    n= len(initial_conditions)
    solution= [[] for i in range(nsteps)]
    solution[0]=[*initial_conditions,ode(0,*initial_conditions)]
    for i in range(1,nsteps):
        mid_solution=euler_method(ode,solution[i-1][:-1],dt/2,1)[1]
        for j in range(1,n+1):
            solution[i].insert(0,taylor(dt,mid_solution[n-j:],j+1))
        solution[i].append(ode(i*dt,*solution[i]))
    return np.array(solution)

def rk4_method(ode,initial_conditions,dt,nsteps):
    nsteps=nsteps+1
    n= len(initial_conditions)
    solution= [[] for i in range(nsteps)]
    solution[0]=[*initial_conditions,ode(0,*initial_conditions)]
    for i in range(1,nsteps):
        prev_t=(i-1)*dt
        prev_x=solution[i-1][0]
        prev_v=solution[i-1][1]
        K1=prev_v
        F1=ode(prev_t,prev_x,K1)
        K2=prev_v+F1*dt/2
        F2=ode(prev_t+dt/2,prev_x+dt*K1/2,K2)
        K3=prev_v+F2*dt/2
        F3=ode(prev_t+dt/2,prev_x+dt*K2/2,K3)
        K4=prev_v+F3*dt
        F4=ode(prev_t+dt,prev_x+dt*K3,K4)
        new_t=prev_t+dt
        new_x=prev_x+(dt/6)*(K1+2*K2+2*K3+K4)
        new_v=prev_v+(dt/6)*(F1+2*F2+2*F3+F4)
        new_a=ode(new_t,new_x,new_v)
        solution[i]=[new_x,new_v,new_a]
    return np.array(solution)

def harmonic_oscilator(_,x0,k=3):
    return -k*x0
def q3_harmonic_oscilator(method):
    ode = lambda t,x,v: harmonic_oscilator(t,x,1)
    nsteps=int(1e5)
    dt=1e-2
    if (method=='euler'):
        nsteps=3968
        dt=2*np.pi/nsteps 
        solution =np.array(euler_method(ode,[1,0],dt,nsteps))
    if (method=='super_euler'):
        nsteps=3968
        dt=2*np.pi/nsteps 
        solution =np.array(super_euler_method(ode,[1,0],dt,nsteps))
    if (method=='rk2'):
        nsteps=34
        dt=2*np.pi/nsteps 
        solution =np.array(rk2_method(ode,[1,0],dt,nsteps))
    if (method=='rk4'):
        nsteps=10
        dt=2*np.pi/nsteps 
        solution =np.array(rk4_method(ode,[1,0],dt,nsteps))
    t_arr = np.array([i*dt for i in range(nsteps+1)])
    if (method=='analytical'):
        plt.plot(t_arr,np.cos(t_arr),label='analytical solution')
        return
    plt.plot(t_arr,solution.transpose()[0],marker='*',ls='--',label='numerical solution ' +method)

def main():
    plt.figure(dpi=1200)
    q3_harmonic_oscilator('euler')
    q3_harmonic_oscilator('super_euler')
    q3_harmonic_oscilator('rk2')
    q3_harmonic_oscilator('rk4')
    q3_harmonic_oscilator('analytical')
    plt.legend()
    plt.grid()

main()

def q3_plot_diff(method):
    ode = lambda t,x,v: harmonic_oscilator(t,x,40)
    nsteps=int(1e5)
    dt=1e-7
    t_arr = np.array([i*dt for i in range(nsteps+1)])
    if (method=='rk4'):
        solution=np.array(rk4_method(ode,[1,0],dt,nsteps))
    analytical=np.cos(np.sqrt(40)*t_arr)
    plt.plot(t_arr,analytical-solution.transpose()[0],label='difference between analytical and numerical')
    plt.grid()
    plt.legend()


def free_fall(t,x,v):
    return -9.8
def q3_free_fall(method):
    ode = free_fall
    nsteps=int(1e4)
    dt=1e-4
    t_arr = np.array([i*dt for i in range(nsteps+1)])
    analytical= 1-0.5*9.8*t_arr**2
    if (method=='euler'):
        solution =np.array(euler_method(ode,[1,0],dt,nsteps))
    if (method=='super_euler'):
        solution =np.array(super_euler_method(ode,[1,0],dt,nsteps))
    if (method=='rk2'):
        solution =np.array(rk2_method(ode,[1,0],dt,nsteps))
    if (method=='rk4'):
        solution =np.array(rk4_method(ode,[1,0],dt,nsteps))
    if (method=='analytical'):
        plt.plot(t_arr,analytical,label='analytical solution')
        return
    plt.plot(t_arr,solution.transpose()[0],label='numerical solution ' +method)

'''
q3_free_fall('euler')
q3_free_fall('super_euler')
q3_free_fall('rk2')
q3_free_fall('analytical')
plt.legend()
plt.grid()

'''

'''
class ode:
    def __init__(self,f,order):
        self.order= order
        self.f = f
'''
