import math
import numpy as np
from sympy import *
import matplotlib.pyplot as plt

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

def optical_bloch_equation(X,gamma=1,omega=1,delta=1):
    rho_gg_dot=(1j*omega/2)*(X[2]-X[3])+gamma*X[1]
    rho_ee_dot=-(1j*omega/2)*(X[2]-X[3])-gamma*X[1]
    rho_ge_dot=-1j*delta*X[2]-(1j*omega/2)*(X[1]-X[0])-(gamma/2)*X[2]
    rho_eg_dot=1j*delta*X[3]+(1j*omega/2)*(X[1]-X[0])-(gamma/2)*X[3]
    return np.array([rho_gg_dot,rho_ee_dot,rho_ge_dot,rho_eg_dot])

def graph_with_name(x,y,name,method,delta_factor):
    plt.figure(dpi=300)
    plt.plot(x,y.real,label=f'real part of {name}')
    plt.plot(x,y.imag,label=f'imaginary part of {name}')
    plt.xlabel('t')
    plt.ylabel(name)
    plt.title(fr'Numerical solution using {method} for $\Delta$={delta_factor}$\cdot$$\Gamma$')
    plt.legend()
    plt.grid()
    plt.show()

def ope_graph(method,gamma=1,delta_factor=1):
    X=[1,0,0,0]
    ode = lambda t,x : optical_bloch_equation(x,gamma=gamma,omega=gamma,delta=gamma*delta_factor)
    nsteps=10000
    dt=10/nsteps
    if (method=='euler'):
        solution =np.array(euler_method(ode,X,dt,nsteps))
    if (method=='rk4'):
        solution =np.array(rk4_method(ode,X,dt,nsteps))
    t_arr = np.array([i*dt for i in range(nsteps+1)])
    graph_with_name(t_arr,solution.transpose()[0],r'$\rho_{gg}$',method,delta_factor)
    graph_with_name(t_arr,solution.transpose()[1],r'$\rho_{ee}$',method,delta_factor)
    graph_with_name(t_arr,solution.transpose()[2],r'$\rho_{ge}$',method,delta_factor)
    graph_with_name(t_arr,solution.transpose()[3],r'$\rho_{eg}$',method,delta_factor)
    #varify_trace
    plt.figure(dpi=300)
    plt.plot(t_arr,solution.transpose()[0].real+solution.transpose()[1].real,label='real part')
    plt.plot(t_arr,solution.transpose()[0].imag+solution.transpose()[1].imag,label='imaginary part')
    plt.title(r'$\rho_{gg}$+$\rho_{ee}$')
    plt.xlabel('t')
    plt.ylabel(r'$\rho_{gg}+\rho_{ee}$')
    plt.legend()
    plt.grid()
    plt.show()

    print(f'In the last iteration rho is:{solution[-1]}')
    print(solution[-1])

def main():
    gamma=3
    ope_graph('rk4',gamma=gamma,delta_factor=1)
    ope_graph('rk4',gamma=gamma,delta_factor=10)
    #ope_graph('euler',gamma=gamma,delta_factor=1)
    #ope_graph('euler',gamma=gamma,delta_factor=10)
    
main()

# solution for steady-state
def steady_state_solution():
    import sympy as sp
    w, x, y, z = sp.symbols('w,x, y, z')
    a1, a2, a3 = sp.symbols('a1:4')
    eq1 = sp.Eq(a1*w - (1j*a2/2)*x + (1j*a2/2)*y , 0)        
    eq2 = sp.Eq((-1j*a2/2)*w-(a1/2+1j*a3)*y+ (1j*a2/2)*z, 0)
    eq3 = sp.Eq((1j*a2/2)*w+(-a1/2+1j*a3)*x+(-1j*a2/2)*z,0)
    eq4 = sp.Eq(w+z,1)
    ans = sp.solve((eq1, eq2, eq3, eq4), (w,x, y, z))
    print(ans)