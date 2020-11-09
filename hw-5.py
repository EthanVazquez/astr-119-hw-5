#!/usr/bin/env python
# coding: utf-8

# In[ ]:


get_ipython().run_line_magic('matplotlib', 'inline')
import numpy as np
import matplotlib.pyplot as plt


# In[ ]:


def func(x):
    a = 2
    b = 10 
    return np.exp(-a*x)*np.cos(b*x) 


# In[ ]:


def func_intergral(x):
    a = 2
    b = 10
    return (5*np.exp(-a*x)*np.sin(b*x))/52. - np.exp(-a*x)*np.cos(b*x)/52. 


# In[ ]:


def trapezoid_core(f,x,h):
    return 0.5*h*( f(x+h) + f(x))


# In[ ]:


def trapezoid_method(f,a,b,N):
    
    x = np.linspace(a,b, N)
    h = x[1]-x[0]
    
    Fint = 0.0
    
    for i in range(0,len(x)-1,1):
        Fint += trapezoid_core(f,x[i],h)
        
    return Fint


# In[ ]:


def simpson_core(f,x,h):
    return h*( f(x) + 4*f(x+h) + f(x+2*h))/3.


# In[ ]:


def simpson_method(f,a,b,N):
    x = np.linspace(a,b,N)
    h = x[1]-x[0]
    
    Fint = 0.0
    
    for i in range(0,len(x)-2,2):
        Fint += simpson_core(f,x[i],h)
        
        
    if((N%2)==0):
        Fint += simpson_core(f,x[-2],0.5*h)
        
    return Fint


# In[ ]:


def romberg_core(f,a,b,i):
    
    h = b-a
    dh = h/2.**(i)
    
    k = h/2.**(i+1)
    
    M = 0.0
    for j in range(2**i):
        M += f(a + 0.5*dh + j*dh)
        
    return k*M


# In[ ]:


def romberg_intergration(f,a,b,tol):
    
    i = 0
    
    imax = 1000
    
    delta = 100.0*np.fabs(tol)
    
    I = np.zeros(imax,dtype=float)
    
    I[0] = 0.5*(b-a)*(f(a) + f(b))
    
    i += 1
    
    while(delta>tol):
        
        I[i] = 0.5*I[i-1] + romberg_core(f,a,b,i)
        
        delta = np.fabs( (I[i]-I[i-1])/I[1] )
        
        print(i,I[i],I[i-1],delta)
        
        if(delta>tol):
            
            i+=1
            
            if(i>imax):
                print("Max iterations reached.")
                raise StopIteration('Stopping iterations after',i)
                
    return I[i]


# In[ ]:


Answer = func_intergral(1)-func_intergral(0)
print(Answer)
print("Trapezoid")
print(trapezoid_method(func,0,1,10))
print("Simpson's Method")
print(simpson_method(func,0,1,10))
print("Romberg")
tolerance = 1.0e-6
RI = romberg_intergration(func,0,1,tolerance)
print(RI, (RI-Answer)/Answer, tolerance)


# In[ ]:




