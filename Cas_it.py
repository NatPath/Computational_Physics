import math
import numpy as np
from sympy import *
import matplotlib.pyplot as plt

def Q3_hell():
    x,v,w,t=symbols('x v w t')
    E_prime=w**2*(x+t*(v-t*w**2*x/2)**2)+(v-t*w**2*(x+t*v/2))**2/2
    E=(w**2*x**2+v**2)/2
    return simplify(E_prime/E)

Q3_hell()