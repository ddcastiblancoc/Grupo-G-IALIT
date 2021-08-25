# -*- coding: utf-8 -*-
"""
Created on Sun Jul 11 14:38:46 2021

@author: Estandar
"""

import numpy as np
import scipy.special as special
import scipy.integrate as integrate
from scipy.optimize import fsolve
from scipy.misc import derivative
from scipy.interpolate import interp1d

#definición de la función de distribución estacionaria
def maxwell_juttner(u,theta):
    return np.heaviside(u,0)*(1/(theta*special.kv(2,1/theta)))*np.exp(-(np.sqrt(1+u**2))/theta)*u**2

#definición de la función de distribución desplazada normalizada numéricamente
def maxwell_juttner_shifted(u,theta,beta):
    gamma = 1/np.sqrt(1-beta**2)
    non_normalized = lambda x: np.heaviside(x,0)*np.exp(-gamma*(np.sqrt(1+x**2))/theta)*x*np.sinh(gamma*beta*x/theta)
    normalization = integrate.quad(non_normalized,0,100*theta)[0]
    return non_normalized(u)/normalization

#----------------------------------------------------------------------------------------------------------------------
#muestreo de la variable aleatoria u en el caso estacionario por el método de rechazo para funciones log-cóncavas
#----------------------------------------------------------------------------------------------------------------------

#la función recibe una temperatura y la cantidad (máxima) de números aleatorios a generar; retorna un array de velocidades
#distribuidas por Maxwell-Juttner
def rejection_log_maxwell_juttner(theta,n):
    
    #moda de la distribución de Maxwell-Juttner estacionaria
    u_m = np.sqrt(2*theta*(theta+np.sqrt(1+theta**2)))

    #resolución numérica de las coordenadas x de las intersecciones de f(u_m)/e con la distribución
    u_minus,u_plus = fsolve(lambda x: np.log(maxwell_juttner(x,theta)/maxwell_juttner(u_m,theta))+1,[u_m/2,3*u_m/2])

    #definición de los lambda y los q
    lambda_plus = -maxwell_juttner(u_plus,theta)/derivative(lambda x: maxwell_juttner(x,theta),u_plus,dx=1e-6)
    lambda_minus = maxwell_juttner(u_minus,theta)/derivative(lambda x: maxwell_juttner(x,theta),u_minus,dx=1e-6)
    q_plus = lambda_plus/(u_plus-u_minus)
    q_minus = lambda_minus/(u_plus-u_minus)
    q_m = 1-(q_plus+q_minus)

    #variables aleatorias uniformes
    U = np.random.uniform(0,1,n)
    V = np.random.uniform(0,1,n)
    
    #inicialización de X
    X = []

    #método de rechazo al muestrear valores de X sobre la función de prueba de colas exponenciales
    for i in range(n):
        #el primer if muestrea alrededor del máximo
        if U[i]<=q_m:
            y = U[i]/q_m
            x = (1-y)*(u_minus+lambda_minus)+y*(u_plus-lambda_plus)
            #estos if aceptam los valores de X que superen la condición de rechazo
            if V[i]<=(maxwell_juttner(x,theta)/maxwell_juttner(u_m,theta)):
                X.append(x)
                continue
        #el segundo if muestrea en la cola exponencial más positiva
        elif U[i]<=(q_m+q_plus):
            e = -np.log((U[i]-q_m)/q_plus)
            x = u_plus-lambda_plus*(1-e)
            if V[i]<=(np.exp(e)*maxwell_juttner(x,theta)/maxwell_juttner(u_m,theta)):
                X.append(x)
                continue
        #el else muestrea en la cola exponencial más cercana al cero
        else:
            e = -np.log((U[i]-(q_m+q_plus))/q_minus)
            x = u_minus+lambda_minus*(1-e)
            if V[i]<=(np.exp(e)*maxwell_juttner(x,theta)/maxwell_juttner(u_m,theta)):
                X.append(x)
                continue          

    #resultado final del muestreo de u
    return np.array(X)

#----------------------------------------------------------------------------------------------------------------------
#muestreo de la variable aleatoria u en el caso estacionario por el método de rechazo para funciones log-cóncavas
#----------------------------------------------------------------------------------------------------------------------

#la función recibe una temperatura y la cantidad (máxima) de números aleatorios a generar; retorna un array de velocidades
#distribuidas por Maxwell-Juttner
def sobol_maxwell_juttner(theta,n):
    
    #variables aleatorias uniformes 1-4
    X1 = np.random.uniform(0,1,n)
    X2 = np.random.uniform(0,1,n)
    X3 = np.random.uniform(0,1,n)
    X4 = np.random.uniform(0,1,n)

    #definición de u y eta según el algoritmo de Sobol
    u = -theta*np.log(X1*X2*X3)
    eta = -theta*np.log(X1*X2*X3*X4)

    #filtro de los valores de u (para la distribución estacionaria) que pasan la condición de Sobol
    u = np.array([u[i] for i in range(0,n) if eta[i]**2-u[i]**2>1])

    return u

#----------------------------------------------------------------------------------------------------------------------
#muestreo de la variable aleatoria u en el caso estacionario por el algoritmo de la transformada inversa
#----------------------------------------------------------------------------------------------------------------------

#la función recibe una temperatura y la cantidad de números aleatorios a generar; retorna un array de velocidades distribuidas
#por Maxwell-Juttner
def inverse_maxwell_juttner(theta,n):
    
    #función de distribución acumulada F de la distribución de interés (se halla integrando numéricamente, ya está normalizada)
    F_maxwell_juttner = lambda u,theta: integrate.quad(lambda x: maxwell_juttner(x,theta),0,u)[0]

    #variable aleatoria uniforme
    X3 = np.random.uniform(0,1,n)

    #rango de u sobre el que hallar los valores de F(u) (lo suficientemente grande para que la interpolación no deje por fuera 
    #muchos valores de X3 al evaluarla)
    u = np.linspace(0,20*theta,2000)
    #valores de F(u)
    F = np.array([F_maxwell_juttner(i,theta) for i in u])

    #interpolación de los valores de la inversa de F al evaluarla en X3 (los valores de X3 que queden por fuera son extrapolados)
    u = interp1d(F,u,fill_value="extrapolate")(X3)
    
    return u

#----------------------------------------------------------------------------------------------------------------------
#muestreo de la magnitud y las componentes de la velocidad por el flipping method para el caso desplazado
#----------------------------------------------------------------------------------------------------------------------

#la función recibe una muestra de velocidades u en el marco en reposo, temperatura, la velocidad del boost al marco de referencia 
#desplazado y la cantidad (máxima) de números aleatorios a generar; retorna los array de las componentes y magnitud de la velocidades 
#distribuida por el Maxwelliano desplazado
def flipping_maxwell_juttner_shifted(u,beta,n):
    
    #factor de Lorentz asociado al boost
    gamma=1/np.sqrt(1-beta**2)

    #variables aleatorias uniformes 5-7 (la cantidad de números aleatorios es la misma de los u estacionarios)
    X5 = np.random.uniform(0,1,u.size)
    X6 = np.random.uniform(0,1,u.size)
    X7 = np.random.uniform(0,1,u.size)

    #componentes de u proyectándolas sobre una esfera de radio u con ayuda de dos variables aleatorias uniformes
    u_x = u*(2*X5-np.full(u.size,1))
    u_y = 2*u*np.sqrt(X5*(np.full(u.size,1)-X5))*np.cos(2*np.pi*X6)
    u_z = 2*u*np.sqrt(X5*(np.full(u.size,1)-X5))*np.sin(2*np.pi*X6)

    #flipping method (los u_x que cumplan una condición cambian de signo, si no, se quedan igual)
    u_x = np.array([-u_x[i] if (-beta*u_x[i]/np.sqrt(1+u[i]**2))>X7[i] else u_x[i] for i in range(0,u.size)])
    #boost en dirección -x (u_x se reasigna, u_y y u_z quedan igual por el boost)
    u_x = np.array([gamma*(u_x[i]+beta*np.sqrt(1+u[i]**2)) for i in range(0,u.size)])

    #nueva magnitud de la velocidad por el boost
    u = np.sqrt(u_x**2+u_y**2+u_z**2)

    return u_x,u_y,u_z,u