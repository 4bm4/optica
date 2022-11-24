from sympy import *
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from numpy.lib.function_base import append
from mpl_toolkits.mplot3d import Axes3D
from sympy.abc import x, y, z
from mpl_toolkits.axisartist.parasite_axes import HostAxes, ParasiteAxes


x, y, z = symbols('x y z')
i, j, k = symbols('i j k')

def conv_fun (fun,var):
    f = lambdify(var, fun)
    return(f)

def parametrizar (vector, punto, min=-0.25,max=0.25):
    if (vector!=None):
        ve1=[]
        ve2=[]
        ve3=[]
        for t in np.arange (min,max,0.5):
            ve1.append(vector[0]*t+punto[0])
            ve2.append(vector[1]*t+punto[1])
            ve3.append(vector[2]*t+punto[2])
        
        return(ve1,ve2,ve3)
    else:
        return(False)

def snell (a,b,c,mu=1.1):
    from sympy.abc import q,w,e

    nueva_matriz=Matrix(( (0, -c,b,mu*b) , (c,0,-a,-a*mu ), (-b,a,0,0) ))
    sist_ec=solve_linear_system(nueva_matriz,q,w,e)

    if (sist_ec!=None):
        if (len(sist_ec)>2):
            vect_final=[sist_ec[q],sist_ec[w],sist_ec[e]]        
            return(vect_final)

def vector_normal(fun, var,a,b):
    
    det=(1+diff(fun,x)**2+diff(fun,y)**2)**(1/2)
    
    vect_nor_x=conv_fun(-diff(fun,x)/det,var)
    vect_nor_y=conv_fun(-diff(fun,y)/det,var)
    vect_nor_z=conv_fun(1/det,var)
    vectfinal=[vect_nor_x(a,b),vect_nor_y(a,b),vect_nor_z(a,b)]
 
    return(vectfinal)

def plot_vect_norm(fun, var,long=0.15, min=-1,max=1,int=0.1):
    
    ax = plt.figure(figsize=(10, 10)).add_subplot(projection='3d')
    
    a, b = np.meshgrid(np.arange(min,max, int),
                      np.arange(min ,max, int))
   
    f_inicial=conv_fun(fun,var)
    
    i,j,k=vector_normal(fun, var,a,b)
    
    #posicion
    puntos_fun_mesh=f_inicial(a,b)
    
    #vecotor normal en función
  
    ax.plot_surface(a,b,puntos_fun_mesh,color='r')
    
    ax.quiver(a, b,puntos_fun_mesh, i, j, k, length=long)
    
    plt.show()
    
    return()



def lineas_normales(fun,var,min=-5,max=5,int=1,alt1=30):
  
    ax = plt.figure(figsize=(10, 10)).add_subplot(projection='3d')
    

    a, b = np.meshgrid(np.arange(min,max, int),
                      np.arange(min ,max, int))
    
    f_inicial=conv_fun(fun,var)
    
    puntos_fun_z=f_inicial(a,b)
    
    ax.plot_surface(a,b,puntos_fun_z,color='r')
    
        
      
    for A in range(-5,5):
        for B in range(-5,5):
            z=f_inicial(A,B)
            z1=[z,z+alt1]
            lin_x=[A,A]
            lin_y=[B,B]
            ax.plot3D(lin_x,lin_y,z1,color='green')

    for A in np.arange (-1,1,0.1):
        for B in np.arange (-1,1,0.1):
            aux=vector_normal(fun, var,A,B)
            vect_salida= snell (aux[0],aux[1],aux[2]) 
            if (parametrizar(vect_salida,[A,B,f_inicial(A,B)])!=False):
                lin_x2,lin_y2,lin_z2=parametrizar(vect_salida,[A,B,f_inicial(A,B)])
                ax.plot3D(lin_x2,lin_y2,lin_z2,color='blue')
    

    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
   

    plt.show()
    
    return()