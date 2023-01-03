import math
from matplotlib import projections
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
import time
import timeit

def c_s(rho):
    return np.sqrt((5/3)*(rho**2/3))

def expansion_of_gas_in_container_initial_conditions(x):
    rho=0
    v=0
    if 0.9<=x<=1:
        rho=1
    if 0<=x<0.9:
        rho=0
    return rho,v 

def calculate_next_rho(lr,v,Nx,dx,dt,order,boundary,g=-10):

    next_rho=[]
    if order=='first':
        if boundary=='periodic':
            for i in range(Nx):
                ip=(i+1)%Nx
                im=(i-1)%Nx
                next_rho_val_i=max(0,(lr[ip]+lr[im])/2-(dt/(2*dx))*(lr[ip]*v[ip]-lr[im]*v[im]))
                #next_rho_val_i=max(0,lr[i]-(dt/(2*dx))*(lr[(i+1)%Nx]*v[(i+1)%Nx]-lr[(i-1)%Nx]*v[(i-1)%Nx]))
                next_rho.append(next_rho_val_i)
        if boundary=='walls':
            for i in range(Nx):
                
                if i==0:
                    next_rho_val_i=max(0,(lr[i]+lr[i+1])/2-(dt/(2*dx))*(lr[i+1]*v[i+1]+lr[i]*v[i]))
                elif i==(Nx-1):
                    next_rho_val_i=max(0,(lr[i]+lr[i-1])/2+(dt/(2*dx))*(lr[i]*v[i]+lr[i-1]*v[i-1]))
                else:
                    next_rho_val_i=max(0,(lr[i+1]+lr[i-1])/2-(dt/(2*dx))*(lr[i+1]*v[i+1]-lr[i-1]*v[i-1]))
                next_rho.append(next_rho_val_i)
    if order=='second':
        if boundary=='periodic':
            for i in range(Nx):
                ip=(i+1)%Nx
                im=(i-1)%Nx
                a=5/3
                zero_order_term= (lr[ip]+lr[im])/2
                first_order_term= -(dt/(2*dx))*(lr[ip]*v[ip]-lr[im]*v[im])
                second_order_term_1=((dt**2)/(2*(dx**2)))*(lr[ip]*v[ip]**2+lr[ip]**a+lr[im]*v[im]**2+lr[im]**a-2*(lr[i]*v[i]**2+lr[i]**a))
                second_order_term_2=-((g*(dt**2))/(4*dx))*(lr[ip]-lr[im])
                second_order_term=second_order_term_1+second_order_term_2
                new_rho_val=zero_order_term+first_order_term+second_order_term
                next_rho_val_i=max(0,new_rho_val)
                next_rho.append(next_rho_val_i)
        if boundary=='walls':
            for i in range(Nx):
                ip=(i+1)%Nx
                im=(i-1)%Nx
                if i==0:
                    temp_lr_im=lr[im]
                    lr[im]=0
                    temp_v_im=v[im]
                    v[im]=0

                if i==(Nx-1):
                    temp_lr_ip=lr[ip]
                    lr[ip]=0
                    temp_v_ip=v[ip]
                    v[ip]=0

                a=5/3
                zero_order_term= (lr[ip]+lr[im])/2
                zero_order_term=lr[i]
                if i==0:
                    first_order_term= -(dt/(2*dx))*(lr[ip]*v[ip]+lr[i]*v[i])
                    second_order_term_1=((dt**2)/(2*(dx**2)))*(lr[ip]*v[ip]**2+lr[ip]**a-(lr[i]*v[i]**2+lr[i]**a))
                    second_order_term_2=-((g*(dt**2))/(4*dx))*(lr[ip]+lr[i])
                                                        
                elif i==(Nx-1):
                    first_order_term= (dt/(2*dx))*(lr[im]*v[im]+lr[i]*v[i])
                    second_order_term_1=((dt**2)/(2*(dx**2)))*(lr[im]*v[im]**2+lr[im]**a-(lr[i]*v[i]**2+lr[i]**a))
                    second_order_term_2=((g*(dt**2))/(4*dx))*(lr[i]+lr[im])

                else:
                    first_order_term= -(dt/(2*dx))*(lr[ip]*v[ip]-lr[im]*v[im])
                    second_order_term_1=((dt**2)/(2*(dx**2)))*(lr[ip]*v[ip]**2+lr[ip]**a+lr[im]*v[im]**2+lr[im]**a-2*(lr[i]*v[i]**2+lr[i]**a))
                    second_order_term_2=-((g*(dt**2))/(4*dx))*(lr[ip]-lr[im])
                second_order_term=second_order_term_1+second_order_term_2
                new_rho_val=zero_order_term+first_order_term+second_order_term
                next_rho_val_i=max(0,new_rho_val)
                next_rho.append(next_rho_val_i)

                if i==0:
                    lr[im]=temp_lr_im
                    v[im]=temp_v_im
                if i==(Nx-1):
                    lr[ip]=temp_lr_ip
                    v[ip]=temp_v_ip

    return next_rho

def calculate_next_v(lr,nr,v,Nx,dx,dt,order,boundary,g=-10,cs=np.sqrt(5/3)):
    next_v=[]

    if boundary=='periodic':
        if order=='first':
            for i in range(Nx):
                a=5/3
                ipp=(i+2)%Nx
                ip=(i+1)%Nx
                im=(i-1)%Nx
                imm=(i-2)%Nx
                if nr[i]==0:
                    next_v.append(0)
                else:
                    new_v=(lr[i]*v[i]-(dt/(2*dx))*(lr[ip]*(v[ip]**2)+lr[ip]**a-(lr[im]*(v[im]**2)+lr[im]**a))+dt*g*lr[i])/nr[i]
                    if np.abs(new_v)>=cs:
                        new_v=np.sign(new_v)*(np.abs(v[i])+cs)/2
                    next_v.append(new_v)

                    #next_v.append((lr[i]*v[i]-(dt/(2*dx))*(lr[ip]*(v[ip]**2)+lr[ip]**(5/3)-(lr[im]*(v[im]**2)+lr[im]**(5/3)))+dt*g*lr[i])/nr[i])
        if order=='second':
            for i in range(Nx):
                # for ease of debugging, i'll use temporary variables
                ipp=(i+2)%Nx
                ip=(i+1)%Nx
                im=(i-1)%Nx
                imm=(i-2)%Nx

                if nr[i]==0:
                    next_v.append(0)
                else:
                    zero_order_term=lr[i]*v[i]
                    first_order_term=(-(dt/(2*dx))*(lr[ip]*(v[ip]**2)+lr[ip]**(5/3)-(lr[im]*(v[im]**2)+lr[im]**(5/3)))+dt*g*0.5*(lr[ip]+lr[im]))
                    a=5/3
                    A1=v[ip]*(lr[ipp]*v[ipp]**2+lr[ipp]**a-(lr[i]*v[i]**2+lr[i]**a))
                    A2=v[im]*(lr[i]*v[i]**2+lr[i]**a-(lr[imm]*v[imm]**2+lr[imm]**a))
                    A=A1-A2
                    B1=(a*lr[ip]**(a-1)-v[ip]**2)*((lr[ipp]*v[ipp]-lr[i]*v[i]))
                    B2=(a*lr[im]**(a-1)-v[im]**2)*((lr[i]*v[i]-lr[imm]*v[imm]))
                    B=B1-B2
                    second_order_term1=((dt**2)/(8*(dx**2)))*(A+B)
                    second_order_term2=-((dt**2)*g/(2*dx))*(lr[ip]*v[ip]-lr[im]*v[im])
                    second_order_term=second_order_term1+second_order_term2
                    new_v=(zero_order_term+first_order_term+second_order_term)/nr[i]
                    if np.abs(new_v)>=cs:
                        new_v=np.sign(new_v)*(np.abs(v[i])+cs)/2
                    next_v.append(new_v)

    if boundary=='walls':
        if order=='first':
            for i in range(Nx):
                a=5/3
                ip=(i+1)%Nx
                im=(i-1)%Nx
                if i==0:
                    temp_lr_im=lr[im]
                    lr[im]=0
                    temp_v_im=v[im]
                    v[im]=0

                if i==(Nx-1):
                    temp_lr_ip=lr[ip]
                    lr[ip]=0
                    temp_v_ip=v[ip]
                    v[ip]=0

                if nr[i]==0:
                    next_v.append(0)
                else:
                    new_v=(lr[i]*v[i]-(dt/(2*dx))*(lr[ip]*(v[ip]**2)+lr[ip]**a-(lr[im]*(v[im]**2)+lr[im]**a))+dt*g*lr[i])/nr[i]
                    if np.abs(new_v)>=cs:
                        new_v=np.sign(new_v)*(np.abs(v[i])+cs)/2
                    next_v.append(new_v)

                if i==0:
                    lr[im]=temp_lr_im
                    v[im]=temp_v_im
                if i==(Nx-1):
                    lr[ip]=temp_lr_ip
                    v[ip]=temp_v_ip

        if order=='second':
            for i in range(Nx):
                ipp=(i+2)%Nx
                ip=(i+1)%Nx
                im=(i-1)%Nx
                imm=(i-2)%Nx
                if i<=1:
                    temp_lr_imm=lr[imm]
                    lr[imm]=0
                    temp_v_imm=v[imm]
                    v[imm]=0
                    if i==0:
                        temp_lr_im=lr[im]
                        lr[im]=0
                        temp_v_im=v[im]
                        v[im]=0

                if i>=(Nx-2):
                    temp_lr_ipp=lr[ipp]
                    lr[ipp]=0
                    temp_v_ipp=v[ipp]
                    v[ipp]=0
                    if i==(Nx-1):
                        temp_lr_ip=lr[ip]
                        lr[ip]=0
                        temp_v_ip=v[ip]
                        v[ip]=0

                if nr[i]==0:
                    next_v.append(0)
                else:
                    '''
                    #zero_order_term=(lr[ip]*v[ip]+lr[im]*lr[im])/(2*nr[i])
                    zero_order_term=(lr[i]*v[i])/nr[i]
                    first_order_term=-(dt/(2*dx))*(lr[ip]*(v[ip]**2)+lr[ip]**(5/3)-(lr[im]*(v[im]**2)+lr[im]**(5/3)))+dt*g*0.5*(lr[i]+nr[i])/nr[i]
                    nd_order_1=(v[ip]**2-(5/3)*lr[ip]**(2/3))*(lr[ipp]*v[ipp]-lr[i]*v[i])
                    nd_order_2=-2*v[ip]*(lr[ipp]*v[ipp]**2+lr[ipp]**(5/3))
                    nd_order_3=-(lr[i]*v[i]**2+lr[i]**(5/3))
                    nd_order_1_3=nd_order_1+nd_order_2+nd_order_3
                    nd_order_4=(v[im]**2-(5/3)*lr[im]**(2/3))*(lr[imm]*v[imm]-lr[i]*v[i])
                    nd_order_5=-2*v[im]*(lr[imm]*v[imm]**2+lr[imm]**(5/3))
                    nd_order_6=-(lr[i]*v[i]**2+lr[i]**(5/3))
                    nd_order_4_6=nd_order_4+nd_order_5+nd_order_6
                    second_order_term=-((dt**2)/(8*dx**2))*(nd_order_1_3-nd_order_4_6)/nr[i]
                    new_v=zero_order_term+first_order_term+second_order_term
                    '''
                    #zero_order_term=(lr[ip]*v[i]+lr[im]*v[im])/(2*nr[i])
                    zero_order_term=lr[i]*v[i]
                    first_order_term=(-(dt/(2*dx))*(lr[ip]*(v[ip]**2)+lr[ip]**(5/3)-(lr[im]*(v[im]**2)+lr[im]**(5/3)))+dt*g*0.5*(lr[i]+nr[i]))
                    nd_order_1=(v[ip]**2-(5/3)*lr[ip]**(2/3))*(lr[ipp]*v[ipp]-lr[i]*v[i])
                    nd_order_2=-2*v[ip]*(lr[ipp]*v[ipp]**2+lr[ipp]**(5/3)-(lr[i]*v[i]**2+lr[i]**(5/3)))
                    #nd_order_3=-(lr[i]*v[i]**2+lr[i]**(5/3))
                    nd_order_1_3=nd_order_1+nd_order_2 #+nd_order_3
                    nd_order_4=(v[im]**2-(5/3)*lr[im]**(2/3))*(lr[imm]*v[imm]-lr[i]*v[i])
                    nd_order_5=-2*v[im]*(lr[imm]*v[imm]**2+lr[imm]**(5/3)-(lr[i]*v[i]**2+lr[i]**(5/3)))
                    #nd_order_6=-(lr[i]*v[i]**2+lr[i]**(5/3))
                    nd_order_4_6=nd_order_4+nd_order_5 #+nd_order_6
                    second_order_term=(-((dt**2)/(8*dx**2))*(nd_order_1_3-nd_order_4_6))
                    new_v=(zero_order_term+first_order_term+second_order_term)/nr[i]
                    if np.abs(new_v)>=cs:
                        new_v=np.sign(new_v)*(np.abs(v[i])+cs)/2
                    next_v.append(new_v)
                    #next_v.append((lr[i]*v[i]-(dt/(2*dx))*(lr[(i+1)%Nx]*(v[(i+1)%Nx]**2)+lr[(i+1)%Nx]**(5/3)-(lr[(i-1)%Nx]*(v[(i-1)%Nx]**2)+lr[(i-1)%Nx]**(5/3)))+dt*g*lr[i])/nr[i])

                if i<=1:
                    lr[imm]=temp_lr_imm
                    v[imm]=temp_v_imm
                    if i==0:
                        lr[im]=temp_lr_im
                        v[im]=temp_v_im

                if i>=(Nx-2):
                    lr[ipp]=temp_lr_ipp
                    v[ipp]=temp_v_ipp
                    if i==(Nx-1):
                        lr[ip]=temp_lr_ip
                        v[ip]=temp_v_ip
        

    return next_v


def hydrodynamics_1d(initial_conditions,boundary='periodic',order='first',xi=0,xf=1,Nx=100,Nt=100,dt=0.0001,g=-10):
    rho=[]
    x_range=np.linspace(xi,xf,Nx)
    dx=(xf-xi)/Nx
    rho.append([initial_conditions(x)[0] for x in x_range])
    v=[]
    v.append([initial_conditions(x)[1] for x in x_range])
    for i in range(Nt):
        rho.append(calculate_next_rho(rho[-1],v[-1],Nx,dx,dt,order,boundary,g))
        v.append(calculate_next_v(rho[-2],rho[-1],v[-1],Nx,dx,dt,order,boundary,g))
    return rho, v

def hw4_q3(order):
    Nx=200
    dx=1/Nx
    cs=np.sqrt(5/3)
    dt=0.02*dx/cs
    Nt=int(np.ceil(3/(dt*cs)))
    g=0
    rho,v = hydrodynamics_1d(expansion_of_gas_in_container_initial_conditions,'periodic',order,0,1,Nx,Nt,dt,g)
    rho = np.array(rho)
    v = np.array(v)
    print(f'total time ran was {Nt*dt}')
    print(f'max of rho is {np.max(rho)}')
    plt.figure(dpi=340)
    xmin,xmax,ymin,ymax=(0,1,0,3)
    plt.imshow(rho, cmap=plt.cm.jet,extent=[xmin,xmax,ymax,ymin], vmin=np.min(rho), vmax=np.max(rho))
    plt.xlabel('x')
    plt.ylabel('t - in scale of ts')
    plt.title(r'heatmap of $\rho$(x,t)')
    plt.grid()
    plt.colorbar() 
    #plt.legend()
    plt.figure(dpi=340)
    xmin,xmax,ymin,ymax=(0,1,0,3)
    plt.imshow(v, cmap=plt.cm.jet,extent=[xmin,xmax,ymax,ymin], vmin=np.min(v), vmax=np.max(v))
    plt.xlabel('x')
    plt.ylabel('t- in the scale of ts')
    plt.title('heatmap of v(x,t)')
    plt.grid()
    plt.colorbar() 
    print(f'In {order} order')
    print("Checking if there is mass conservation:")
    print(f'The difference in the masses between start and finish is {np.sum(rho[-1])-np.sum(rho[0])}')

    v=np.array(v)
    rho=np.array(rho)
    print("Checking if there is momentum conservation :")
    print(f'The difference in the momentum between start and finish is {np.sum(v[-1]*rho[-1])-np.sum(v[0]*rho[0])}')

    sum_rho=np.array([ np.sum(rho[i]*dx) for i in range(Nt)])
    sum_v= np.array([ np.sum(v[i]*dx) for i in range(Nt)])
    plt.figure(dpi=340)
    plt.plot(np.array(list(range(Nt)))*dt,sum_rho)
    plt.title(rf'Total mass by time')
    plt.xlabel('Time')
    plt.ylabel('Total Mass')
    plt.grid()
    plt.show()

    plt.figure(dpi=340)
    sum_rho= [ np.sum(rho[i]) for i in range(Nt)]
    plt.plot(np.array(list(range(Nt)))*dt,sum_rho*sum_v)
    plt.title(rf'Total momentum by time')
    plt.xlabel('Time')
    plt.ylabel('Total momentum')
    plt.grid()
    plt.show()

def analytical_stable_hydrostatic(x):
    if x<0.25:
        return (-4*x+1)**(3/2)
    else:
        return 0

def hw4_q4(order):
    Nx=200
    dx=1/Nx
    cs=np.sqrt(5/3)
    dt=0.2*dx/cs
    Nt=int(np.ceil(3/(dt*cs)))
    g=-10
    rho,v = hydrodynamics_1d(expansion_of_gas_in_container_initial_conditions,'walls',order,0,1,Nx,Nt,dt,g)
    rho = np.array(rho)
    v = np.array(v)
    print(f'total time ran was {Nt*dt}')
    print(f'max of rho is {np.max(rho)}')
    plt.figure(dpi=340)
    xmin,xmax,ymin,ymax=(0,1,0,3)
    plt.imshow(rho, cmap=plt.cm.jet,extent=[xmin,xmax,ymax,ymin], vmin=np.min(rho), vmax=np.max(rho))
    plt.xlabel('x')
    plt.ylabel('t - in scale of ts')
    plt.title(r'heatmap of $\rho$(x,t)')
    plt.grid()
    plt.colorbar() 
    #plt.legend()
    plt.figure(dpi=340)
    xmin,xmax,ymin,ymax=(0,1,0,3)
    plt.imshow(v, cmap=plt.cm.jet,extent=[xmin,xmax,ymax,ymin], vmin=np.min(v), vmax=np.max(v))
    plt.xlabel('x')
    plt.ylabel('t- in the scale of ts')
    plt.title('heatmap of v(x,t)')
    plt.grid()
    plt.colorbar() 
    print(f'In {order} order')
    print("Checking if there is mass conservation:")
    print(f'The difference in the masses between start and finish is {np.sum(rho[-1])-np.sum(rho[0])}')

    v=np.array(v)
    rho=np.array(rho)
    print("Checking if there is momentum conservation :")
    print(f'The difference in the momentum between start and finish is {np.sum(v[-1]*rho[-1])-np.sum(v[0]*rho[0])}')

    x_range=np.linspace(0,1,Nx)
    hydrostatic_stable=[analytical_stable_hydrostatic(x) for x in x_range]
    plt.figure(dpi=340)
    plt.plot(x_range,rho[-1],label=r'last timestamp $\rho$(x) numerical')
    plt.plot(x_range,rho[0],label=r'first timestamp $\rho$(x) numerical')
    plt.plot(x_range,rho[math.ceil(0.125*Nt/3)],label=r't=1/8 timestamp $\rho$(x) numerical')
    plt.plot(x_range,rho[math.ceil(0.25*Nt/3)],label=r't=1/4 timestamp $\rho$(x) numerical')
    plt.plot(x_range,rho[math.ceil(0.5*Nt/3)],label=r't=1/2 timestamp $\rho$(x) numerical')
    plt.plot(x_range,rho[math.ceil(0.75*Nt/3)],label=r't=3/4 timestamp $\rho$(x) numerical')
    plt.plot(x_range,rho[math.ceil(Nt/3)],label=r't=1 timestamp $\rho$(x) numerical')
    plt.plot(x_range,hydrostatic_stable,label='analytical solution')
    plt.xlabel('x')
    plt.ylabel(r'$\rho$')
    plt.legend()
    plt.title(rf'$\rho$(x)')
    plt.grid()
def main():
    #hw4_q3('first')
    #hw4_q3('second')

    hw4_q4('first')

    #hw4_q4('second')


main()