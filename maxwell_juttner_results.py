# -*- coding: utf-8 -*-
"""
Created on Fri Jul  9 14:21:48 2021

@author: Estandar
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

#introducción de los algoritmos para muestrear el caso estacionario y para boostear al caso desplazado
import maxwell_juttner_algorithms as mj

#texto TeX en gráficas (des-comentar si tiene un compilador TeX instalado, p. ej., Texmaker)
"""
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": "Computer Modern Roman"})
"""

#número (máximo) de valores aleatorios muestreados de los algoritmos
n=1000000

#gráficas de la distribución muestreada comparada a la teórica para el caso estacionario para 3 temperaturas en el array 
#theta (puede agregar más si desea)

theta=np.array([1,100,10000])

plt.figure()
for i in range(0,theta.size):
    #cálculo del histograma de velocidades en el caso estacionario por medio del algoritmo de Sobol
    H = np.histogram(mj.sobol_maxwell_juttner(theta[i],n),100,density=True)
    bin_centers = 0.5*(H[1][1:]+H[1][:-1])
    #3 primeros subplots
    plt.subplot(2,theta.size,i+1)
    #gráfica de la distribución teórica
    plt.plot(bin_centers,mj.maxwell_juttner(bin_centers,theta[i]),'k-',linewidth=1)
    #gráfica de las counts del histograma
    plt.plot(bin_centers,H[0],'r-',linewidth=1,label=r'$\theta=$ '+str(theta[i]))
    #personalización
    plt.legend()
    plt.ticklabel_format(axis='both',style='sci',scilimits=(-2,2))
    plt.xlabel(r'$u$')
    if i==0:
        plt.ylabel(r'$f(u)$ (estacionario)')

#gráficas de la distribución muestreada comparada a la teórica para el caso desplazado para theta fijo (puede cambiarse 
#si desea) y velocidades del boost en el array beta (puede agregar más si desea, debe ser el mismo número que el array 
#theta de arriba)

theta=1
beta=np.array([0.01,0.5,0.95])
#aplicar el boost por el flipping method requiere la correspondiente muestra en el caso estacionario; se escoge por el 
#algoritmo de Sobol
u = mj.sobol_maxwell_juttner(theta,n)

for i in range(0,beta.size):
    #cálculo del histograma de velocidades en el caso desplazado con el flipping method
    H = np.histogram(mj.flipping_maxwell_juttner_shifted(u,beta[i],n)[3],100,density=True)
    bin_centers = 0.5*(H[1][1:]+H[1][:-1])
    #3 siguientes subplots
    plt.subplot(2,beta.size,i+1+beta.size)
    #gráfica de la distribución teórica
    plt.plot(bin_centers,mj.maxwell_juttner_shifted(bin_centers,theta,beta[i]),'k-',linewidth=1)
    #gráfica de las counts del histograma
    plt.plot(bin_centers,H[0],'r-',linewidth=1,label=r'$\theta=$ '+str(theta)+r', $\beta=$ '+str(beta[i]))
    #personalización
    plt.legend()
    plt.ticklabel_format(axis='both',style='sci',scilimits=(-2,2))
    plt.xlabel(r'$u$')
    if i==0:
        plt.ylabel(r'$f(u)$ (desplazado)')

#gráficas del error entre simulado y teórico como la razón entre ambos casos (para el caso desplazado); se toman un theta 
#y beta fijos que pueden cambiarse

theta=1    
beta=0.95
#muestras de las velociadades en el caso estacionario para los tres algoritmos
u_sobol = mj.sobol_maxwell_juttner(theta,n)
u_rejection = mj.rejection_log_maxwell_juttner(theta,n)
u_inverse = mj.inverse_maxwell_juttner(theta,n)        

plt.figure()

#gráfica para el algoritmo de rechazo para funciones log-cóncavas
plt.subplot(1,3,1)
H = np.histogram(mj.flipping_maxwell_juttner_shifted(u_rejection,beta,n)[3],100,density=True)
bin_centers = 0.5*(H[1][1:]+H[1][:-1])
plt.plot(bin_centers,H[0]/mj.maxwell_juttner_shifted(bin_centers,theta,beta),'r-',linewidth=1,label=r'rejection for log-concave algorithm')
#personalización
plt.legend()
plt.axis([0, 60, 0.5, 2])
plt.yscale('log', basey=2)
plt.gca().yaxis.set_major_formatter(mticker.ScalarFormatter())
plt.xlabel(r'$u$')
plt.ylabel(r'$f(u)$ simulado/$f(u)$ teórico')

#gráfica para el algoritmo de Sobol
plt.subplot(1,3,2)
H = np.histogram(mj.flipping_maxwell_juttner_shifted(u_sobol,beta,n)[3],100,density=True)
bin_centers = 0.5*(H[1][1:]+H[1][:-1])
plt.plot(bin_centers,H[0]/mj.maxwell_juttner_shifted(bin_centers,theta,beta),'r-',linewidth=1,label=r'Sobol algorithm')
#personalización
plt.legend()
plt.axis([0, 60, 0.5, 2])
plt.yscale('log', basey=2)
plt.gca().yaxis.set_major_formatter(mticker.ScalarFormatter())
plt.xlabel(r'$u$')

#gráfica para el algoritmo de la transformada inversa
plt.subplot(1,3,3)
H = np.histogram(mj.flipping_maxwell_juttner_shifted(u_inverse,beta,n)[3],100,density=True)
bin_centers = 0.5*(H[1][1:]+H[1][:-1])
plt.plot(bin_centers,H[0]/mj.maxwell_juttner_shifted(bin_centers,theta,beta),'r-',linewidth=1,label=r'inverse transform algorithm')
#personalización
plt.legend()
plt.axis([0, 60, 0.5, 2])
plt.yscale('log', basey=2)
plt.gca().yaxis.set_major_formatter(mticker.ScalarFormatter())
plt.xlabel(r'$u$')