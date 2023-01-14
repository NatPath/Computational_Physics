import math
from matplotlib import projections
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.cm import ScalarMappable

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

def newton_raphson(f,val,intial_guess,precision):
    '''
    Newton Raphson for differantiable functions
    gets a 
        function - f
        value wanted - val
        intial guess of root - intital_guess
        precision wanted, how close the value should be to val - precision
    returns
        (root,f(root))
    '''
    x=symbols('x')
    guess=intial_guess
    func_guess = N(f(x).subs(x,guess))    
    while abs(func_guess-val)>precision:
        diff_f_guess=N(diff(f(x),x).subs(x,guess))
        mid_guess=abs(func_guess-val)
        if diff_f_guess==0:
            print('foo')
            return None,None
        new_guess= guess-(func_guess-val)/diff_f_guess
        guess=new_guess
        func_guess = N(f(x).subs(x,new_guess))   
    root = guess
    return root ,f(root)


def newton_raphson_old(of,val,intial_guess,precision):
    '''
    Newton Raphson for differantiable functions
    gets a 
        function - f
        value wanted - val
        intial guess of root - intital_guess
        precision wanted, how close the value should be to val - precision
    returns
        (root,f(root))
    '''
    x=symbols('x')
    guess=intial_guess
    f=lambda t: of(t)-val
    func_guess = N(f(x).subs(x,guess))    
    while abs(func_guess)>precision:
        diff_f_guess=N(diff(f(x),x).subs(x,guess))
        mid_guess=abs(func_guess-val)
        if diff_f_guess==0:
            print('foo')
            return None,None
        new_guess= guess-func_guess/diff_f_guess
        guess=new_guess
        func_guess = N(f(x).subs(x,new_guess))   
    root = guess
    return root ,f(root)

def t_by_t_ret(t_ret,path,point):
    res=t_ret+sqrt((path(t_ret)[0]-point[0])**2+(path(t_ret)[1]-point[1])**2)
    return res

def plot_t_minus_tret(t_range,observer,v):
    print(f'[t-t_ret](t) observer at (x,y)={observer} v={v}')
    x=observer[0]
    y=observer[1]
    f= lambda t: (t-v*y-np.sqrt((v*y-t)**2+(x**2+y**2-t**2)*(1-v**2)))/(1-v**2)
    analytical_sol=[t-f(t) for t in t_range]
    path = lambda t : np.array([0,v*t])
    t_by_t_ret_specified = lambda t_ret: t_by_t_ret(t_ret,path,np.array([x,y]))
    newton_sol=[t-newton_raphson(t_by_t_ret_specified,t,t,0.0001)[0] for t in t_range]
    plt.plot(t_range,analytical_sol,label=f'analytical- observer={observer} v={v}')
    plt.plot(t_range,newton_sol,label=f'numerical by newton - observer={observer} v={v}')
    plt.title('t-t_ret as a function of t')
    plt.legend()
    plt.xlabel('t')
    plt.ylabel('t-t_ret')
    plt.grid()

def q1():
    t_range=np.arange(-10,10,0.1)
    observer_1=(0.1,0)
    observer_2=(1,0)
    observer_3=(10,0)
    v_1=0.1
    v_2=0.9

    plt.figure(dpi=160)
    plot_t_minus_tret(t_range,observer_1,v_1)
    plot_t_minus_tret(t_range,observer_2,v_1)
    plot_t_minus_tret(t_range,observer_3,v_1)
    plt.show()

    plt.figure(dpi=160)
    plot_t_minus_tret(t_range,observer_1,v_2)
    plot_t_minus_tret(t_range,observer_2,v_2)
    plot_t_minus_tret(t_range,observer_3,v_2)
    plt.show()

def t_ret(t,path,observer,precision=0.0001):
    x,y=observer
    t_by_t_ret_specified = lambda t_ret: t_by_t_ret(t_ret,path,np.array([x,y]))
    res=newton_raphson(t_by_t_ret_specified,t,t,precision)[0]
    return float(res)


def q2():
    x_range=np.arange(-10,10,1)
    w_1=0.1
    w_2=0.9

    path_1=lambda t : np.array([cos(w_1*t),sin(w_1*t)])
    path_2=lambda t : np.array([cos(w_2*t),sin(w_2*t)])

    sol_1=np.array([[t-t_ret(t,path_1,np.array([x,0])) for t in np.arange(0,4*np.pi/w_1,0.3)] for x in x_range])
    print('finished sol_1')
    sol_2=np.array([[t-t_ret(t,path_2,np.array([x,0])) for t in np.arange(0,4*np.pi/w_2,0.03)] for x in x_range])

    np.save(f'hw6_q{2}_omega1.npy',sol_1)
    np.save(f'hw6_q{2}_omega2.npy',sol_2)

    plt.figure(dpi=340)
    xmin,xmax,ymin,ymax=(0,4*np.pi/w_1,-100,100)
    plt.imshow(sol_1, cmap=plt.cm.jet,extent=[xmin,xmax,ymax,ymin])
    plt.xlabel('t')
    plt.ylabel('x [10^-1 scale]')
    plt.title(r'heatmap of [t-t_ret](x,y=0,t) for $\omega$1')
    plt.grid()
    plt.colorbar() 


    plt.figure(dpi=340)
    xmin,xmax,ymin,ymax=(0,4*np.pi/w_2,-10,10)
    plt.imshow(sol_2, cmap=plt.cm.jet,extent=[xmin,xmax,ymax,ymin])
    plt.xlabel('t')
    plt.ylabel('x')
    plt.title(r'heatmap of [t-t_ret](x,y=0,t) for $\omega$2')
    plt.grid()
    plt.colorbar() 

def phi_orbit(t,observer,omega,precision=0.0001):
    x=observer[0]
    y=observer[1]
    path=lambda t : np.array([cos(omega*t),sin(omega*t)])
    t_by_t_ret_specified = lambda t_ret: t_by_t_ret(t_ret,path,np.array([x,y]))
    t_ret=float(newton_raphson(t_by_t_ret_specified,t,t,precision)[0])
    res = 1/(t-t_ret-omega*(y*np.cos(omega*t_ret)-x*np.sin(omega*t_ret)))
    return res

def q3():
    x_range=np.arange(-100,100,1)
    y_range=np.arange(-100,100,1)

    w_1=0.2
    w_2=0.9
    phi_1_sol_t_0=[[phi_orbit(0,(x,y),w_1) for y in y_range] for x in x_range]
    print('finished sol_1')
    phi_2_sol_t_0=[[phi_orbit(0,(x,y),w_2) for y in y_range] for x in x_range]
    np.save(f'hw6_q3_omega1.npy',phi_1_sol_t_0)
    np.save(f'hw6_q3_omega2.npy',phi_2_sol_t_0)

    plt.figure(dpi=340)
    xmin,xmax,ymin,ymax=(-100,100,-100,100)
    plt.imshow(phi_1_sol_t_0, cmap=plt.cm.jet,extent=[xmin,xmax,ymax,ymin])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(fr'heatmap of $\phi$(x,y,t=0) for $\omega$1={w_1}')
    plt.grid()
    plt.colorbar() 


    plt.figure(dpi=340)
    xmin,xmax,ymin,ymax=(-100,100,-100,100)
    plt.imshow(phi_2_sol_t_0, cmap=plt.cm.jet,extent=[xmin,xmax,ymax,ymin])
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(fr'heatmap of $\phi$(x,y,t=0) for $\omega$1={w_2}')
    plt.grid()
    plt.colorbar() 

def E_orbit(t,observer,w,precision=0.0001):
    x=observer[0]
    y=observer[1]
    path=lambda t : np.array([cos(w*t),sin(w*t),0])
    v= lambda t:w*np.array([-np.sin(w*t),np.cos(w*t),0])
    t_by_t_ret_specified = lambda t_ret: t_by_t_ret(t_ret,path,np.array([x,y]))
    t_ret=float(newton_raphson(t_by_t_ret_specified,t,t,precision)[0])
    v=v(t_ret)
    path_tret=path(t_ret)
    path=np.array([float(path_tret[0]),float(path_tret[1]),float(path_tret[2])])
    observer=np.array([x,y,0])
    eta=observer-path
    eta_norm=np.linalg.norm(eta)
    u=eta/eta_norm-v
    a=-w**2*path
    res=(eta_norm/(np.dot(eta,u))**3)*((1-np.dot(v,v))*u+np.cross(eta,np.cross(u,a)))
    return res

def plot_electric_field(Ex, Ey):
    plt.figure(dpi=340)
    x, y = np.linspace(-100, 100, Ex.shape[0]), np.linspace(-100, 100, Ex.shape[1])
    X, Y = np.meshgrid(x, y)

    norm = np.sqrt(Ex**2 + Ey**2)
    #print(norm)
    Ex = Ex / norm
    Ey = Ey / norm

    # Plot the vector field with color parameter
    plt.quiver(X, Y, Ex, Ey ,norm)

    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid()
    plt.show()
    
def q4():
    x_range=np.arange(-100,100,4)
    y_range=np.arange(-100,100,4)

    w_1=0.2
    w_2=0.9
    E_1_sol_t_0=np.array([[E_orbit(0,(x,y),w_1) for y in y_range] for x in x_range])
    E_1_sol_t_0_transposed=np.transpose(E_1_sol_t_0,(2,0,1))
    print('finished sol_1')
    E_2_sol_t_0=np.array([[E_orbit(0,(x,y),w_2) for y in y_range] for x in x_range])
    E_2_sol_t_0_transposed=np.transpose(E_2_sol_t_0,(2,0,1))
    np.save(f'hw6_q4_E_omega1.npy',E_1_sol_t_0)
    np.save(f'hw6_q4_E_omega2.npy',E_2_sol_t_0)

    plot_electric_field(E_1_sol_t_0_transposed[0],E_1_sol_t_0_transposed[1])
    plot_electric_field(E_2_sol_t_0_transposed[0],E_2_sol_t_0_transposed[1])

def main():
    #q1()
    #q2()
    #q3()
    q4()
main()