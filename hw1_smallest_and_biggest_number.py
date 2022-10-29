import math
import numpy as np
from sympy import *
import matplotlib.pyplot as plt
## 1.1
#find the smallest and biggest number the computer can represent
# we'll do a "binary search" for the smallest number which does not yeild equal 0
def smallest_number():
    high=1
    divider=2
    candidate = high*(1-1/divider)
    while high != candidate:
        if candidate>0:
            high = candidate 
        else:
            divider=divider*2
        candidate = high*(1-1/divider)
    return candidate
#print("The smallest positive number represented in python is: " , smallest_number())

# This is an attempt to find the biggest number python
# This code will run untill the memory of the computer it is running on will run out.
# That is because of the way python represents integer numbers
# The esseence of it is that python represents an integers with an array of digits with unlimited size.
# Note that the is a build in value in python of infinity, float('inf')- which is bigger than any integer
def biggest_number():
    candidate=2
    prev_candidate=1
    multiplier = 2
    count=0
    while candidate != prev_candidate:
        count=count+1
        if count==int(1e6):# the loop will run infinitly, 1million iterations are enough to show the point..
            break
        if candidate == float('inf'):
            multiplier = (1+ multiplier)/2
            candidate = prev_candandiate
        prev_candandiate = candidate
        candidate=candidate*multiplier
    return candidate
#print("The biggest number in python is: ", biggest_number())

## 1.2
#find the smallest epsilon s.t 1+epsilon != 1
def smallest_epsilon_rough():
    epsilon = 1
    prev_epsilon =10
    divider=10
    while epsilon!=prev_epsilon:
        if (1+epsilon)==1:
            divider=divider*10
            epsilon = prev_epsilon
        prev_epsilon=epsilon
        epsilon = epsilon*(1-1/divider)
    return prev_epsilon
#print("Smallest epsilon s.t 1+epsilon != 1 is: ", smallest_epsilon_rough())


## 1.3
def nth_term_coefiecient_of_sqrt_x_plus_one(n):
    res=1
    k=1/2
    for _ in range(n):
        res=res*k
        k= k-1
    return res/math.factorial(n)

def sqrt_of_one_point_one():
    res_exp=np.sqrt(1.1)
    res_taylor=0
    i=0
    while np.abs(res_taylor -res_exp)>3e-16:
        res_taylor=res_taylor+nth_term_coefiecient_of_sqrt_x_plus_one(i)*0.1**i
        i=i+1
    return res_taylor,i
'''
res_taylor,i=sqrt_of_one_point_one()
print("Number of terms needed to reach the computer accuracy is: ", i)
print("The result up to that term was: ", res_taylor)
print("In comparision to the value calculated: ", np.sqrt(1.1))
''' 


    
## 2
def f(x):
    return log(x)**2+x

def R(func,func_div,x_val,dx):
    df= func(x_val+dx)-func(x_val)
    func_div_val=func_div(x_val)
    return ((func_div_val-df/dx)/func_div_val)
def plot_R_ln_dx(func,x_val):
    dx_range1=np.array([np.power(10,x) for x in np.linspace(-20,-1,1000000)])
    dx_range2=np.linspace(1,int(1e6),int(1e6))
    dx_range=np.concatenate((dx_range1,dx_range2))
    x = symbols('x')
    func_diff = lambdify(x,diff(func(x),x))
    func = lambdify(x,func(x))
    R_range= np.array([abs(R(func,func_diff,x_val,dx)) for dx in dx_range])
    dx_min_arg = np.argmin(R_range)
    dx_min= dx_range[dx_min_arg]
    R_min=R_range[dx_min_arg]
    print("dx min is: ", dx_min)
    print("dx min arg is: ", dx_min_arg)
    print("|R| min is: ", R_min)
    plt.figure(dpi=1200)
    plt.plot(np.log(dx_range),R_range)
    plt.grid()
    plt.xlabel('log(dx)')
    plt.ylabel('R')
#plot_R_ln_dx(f,3)

#3
def y(x,a,b):
    return a*x**2+b*x

def euler_plot(steps,a,v0_x,v0_y,T,midpoint=False):
    s=""
    if midpoint:
        s="midpoint"
    x_t = lambda t: v0_x*t
    y_0=0
    y_arr_prime=[v0_y]
    y_arr=[y_0]
    dt=T/steps
    for i in range(1,steps+1):
        y_arr_prime.append(a*dt+y_arr_prime[i-1])
        if (midpoint==True):
            y_arr.append(y_arr[i-1]+dt*(y_arr_prime[i-1]+y_arr_prime[i])/2)
        else:
            y_arr.append(y_arr[i-1]+dt*y_arr_prime[i-1])
    fig=plt.plot(np.array([x_t(t) for t in np.linspace(0,T,steps+1)]),y_arr,label="euler {nsteps} steps {nmidpoint}".format(nsteps=steps,nmidpoint=s))
    return fig

def plot_balistic(T,v0_x,v0_y):
    #analytical
    g=9.8
    x_max=T*v0_x
    x_range = np.linspace(0,x_max,1000)
    a=-g/(2*v0_x**2)
    b=v0_y/v0_x
    y_range=np.array([y(x,a,b) for x in x_range])
    plt.plot(x_range,y_range,label='analytical')
    plt.grid()
    plt.xlabel('x')
    plt.ylabel('y')

    euler_plot(10,-g,v0_x,v0_y,T)
    euler_plot(100,-g,v0_x,v0_y,T)
    #midpoint
    euler_plot(10,-g,v0_x,v0_y,T,True) # much smaller time step
    euler_plot(100,-g,v0_x,v0_y,T,True) # baseline time step
    euler_plot(1000,-g,v0_x,v0_y,T,True) # much larger time step
    euler_plot(1,-g,v0_x,v0_y,T,True) # single time step

    plt.legend()
    plt.show() 
plot_balistic(3,1,10)

