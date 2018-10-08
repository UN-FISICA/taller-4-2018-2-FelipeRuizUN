#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  6 19:57:56 2018

@author: oriol
"""
import math as m 
from scipy import optimize  
def dif(f,x,h):
     return (f(x+h/2)-f(x-h/2))/h
    
    
class Derivada:
    def __init__(self,f,metodo='adelante',dx=0.001):
        self.f=f
        self.metodo=metodo
        self.dx=dx
    def calc(self,x):
        
        if self.metodo=='adelante':
            der=(self.f(x+self.dx)-self.f(x))/self.dx
            return der
        elif self.metodo=='central':
            der=dif(self.f,x,self.dx)
            return der
        elif self.metodo=='extrapolada':
            der=(dif(self.f,x,self.dx/2)-dif(self.f,x,self.dx))/2
            return der            
        elif self.metodo=='segunda':
            der=(self.f(x+self.dx)-self.f(x-self.dx)-2*self.f(x))/(self.dx**2)
            return der
        else:
            return NotImplemented
class Zeros:
    def __init__(self,f,metodo,error=1e-4,max_iter=100):
        self.f=f
        self.metodo=metodo
        self.error=error
        self.max_iter=max_iter 
    def zero(self,vi):
        if self.metodo=='newton':
            x0=vi
            i=0
            a=Derivada(self.f,'extrapolada',1e-10)
            while (self.f(x0)>=self.error) and (i<=self.max_iter):
                x0=x0-self.f(x0)/a.calc(x0)
                i+=1
                return x0
        elif self.metodo=='bisectriz':
            a=vi[0]
            b=vi[1]
            x0=a
            i=0
            while (abs(self.f(x0))>=self.error) and (i<=self.max_iter):
                x0=(a+b)/2
                if self.f(x0)*self.f(a)>0:
                    a=x0
                elif self.f(x0)*self.f(b)>0:
                    b=x0
                elif self.f(x0)==0:
                    break
                i+=1
            return x0         
        elif self.metodo=='interpolacion':
            a=vi[0]
            b=vi[1]
            x0=a
            i=0
            while (abs(self.f(x0))>=self.error) and (i<=self.max_iter):
                x0=a-self.f(a)*(b-a)/(self.f(b)-self.f(a))
                if self.f(x0)*self.f(a)>0:
                    a=x0
                elif self.f(x0)*self.f(b)>0:
                    b=x0
                elif self.f(x0)==0:
                    break
                i+=1
            return x0
        elif self.metodo=='fsolve-sp': 
            return optimize.fsolve(self.f,vi,xtol=self.error,maxfev=self.max_iter) 
        elif self.metodo=='brentq-sp':
            return optimize.brentq(self.f,vi[0],vi[1],xtol=self.error,maxiter=self.max_iter)
        else:
            return NotImplemented
if __name__ == "__main__":
    c=Derivada(m.sin,'central',1e-7)
    print(c.calc(0))
    d=Zeros(m.sin,'interpolacion',1e-7,200)
    print(d.zero([3,3.3]))
    
        
    
            
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
