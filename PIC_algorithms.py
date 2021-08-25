# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 23:24:06 2021

@author: Estandar
"""

import numpy as np
import pandas as pd

#importación del muestreo de u iniciales por la distribución de Maxwell-Juttner por el algoritmo de Sobol
import maxwell_juttner_algorithms as mj

#implementación de la proyección directa para la densidad de carga inicial de las posiciones de las partículas del IGM x_p,y_p 
#(de tamaño de celda 2d) en el punto de grilla x,y
def projection(x,y,x_p,y_p,d):
    #filtros para encontrar las partículas cuyas celdas encierran el punto de grilla dado x,y
    filter_x = np.array([i for i in range(x_p.size) if np.abs(x-x_p.loc[i])<=d])
    filter_y = np.array([i for i in filter_x if np.abs(y-y_p.loc[i])<=d])
    
    #si no hay ninguna partícula que encierre el punto de grilla, retorna cero la densidad de carga; si hay partículas, 
    #suma sobre todas las partículas (la mitad son iones, y la otra mitad electrones), evaluando en el punto de grilla 
    #cada función de forma
    if filter_y.size!=0:
        P=0
        for i in filter_y:
            if i<x_p.shape[0]/2:
                P+=(1-np.abs((x-x_p.loc[i])/d))*(1-np.abs((y-y_p.loc[i])/d))/8
            else:
                P+=-(1-np.abs((x-x_p.loc[i])/d))*(1-np.abs((y-y_p.loc[i])/d))/8
        return P
    else:
        return 0
    
#algoritmo de sobre-relajación para resolver la ecuación de Poisson con la distribución de carga inicial rho sacada de proyección
def overrelaxation(N,h,X,Y,rho):
    tol=1e-5
    omega=1.5
    u2=np.zeros([N,N])
    u=u2.copy()
    
    #se calcula el potencial en cada punto de la grilla X,Y (NxN) utilizando diferencias finitas, pero se le agrega 
    #un término adicional con u2 y 1<omega<2 para más rápida convergencia; el while para cuando el error definido (u con u2) 
    #es más pequeño que la tolerancia
    while True:
        for i in range(N):
            for j in range(N):
                if 0<i<N-1 and j==0:
                    u[i,0]=(1-omega)*u2[i,0]+omega*(u[i-1,0]+u[i+1,0]+u[i,N-2]+u[i,1]+rho.iloc[i,0]*h**2)/4
                    u[i,N-1]=u[i,0]
                if i==0 and 0<j<N-1:
                    u[0,j]=(1-omega)*u2[0,j]+omega*(u[N-2,j]+u[1,j]+u[0,j-1]+u[0,j+1]+rho.iloc[0,j]*h**2)/4
                    u[N-1,j]=u[0,j]
                elif 0<i<N-1 and 0<j<N-1:
                    u[i,j]=(1-omega)*u2[i,j]+omega*(u[i-1,j]+u[i+1,j]+u[i,j-1]+u[i,j+1]+rho.iloc[i,j]*h**2)/4
        err=np.max(np.abs(u-u2))
        if err<tol:
            break
        u2=u.copy()

    #se calcula el campo eléctrico utilizando la función gradient y el paso de grilla h para diferencias finitas; se retornan 
    #los campos como DataFrames
    E_x,E_y=np.gradient(u,h)
    E_x=-E_x
    E_y=-E_y
    
    return pd.DataFrame(E_x),pd.DataFrame(E_y)

#perfil del campo magnético inicial que sólo depende de y con ancho W_b
def B_z_initial(X,Y,W_b):
    X,Y=np.meshgrid(X,Y)    
    return pd.DataFrame(np.tanh(Y/W_b)).T

#el algoritmo de Sobol para cargar velocidades por Maxwell-Juttner tiene cierta probabilidad de rechazar un intento; 
#si se quieren obtener n velocidades, es posible que se obtengan menos; es necesario hacer una repetición del intento hasta que
#se obtengan los n valores exactos
def sobol_filter(theta,n):
    while True:
        u_x_stationary = mj.sobol_maxwell_juttner(theta,n)
        if u_x_stationary.size==n:
            break
    return u_x_stationary

#implementación de la interpolación de una función Z evaluada en los puntos de grilla X,Y, en un punto intermedio x,y 
#(posiciones de las partículas específicamente)
def interpolation(x,y,X,Y,Z,d):
    #filtros que buscan el punto de grilla más cercano al punto x,y, o más específicamente, los índices de los valores máximos 
    #de X,Y tales que X<=x y Y<=y
    for k in range(X.size):
        if X[k]>x:
            i=k-1
            break
        else:
            i=k
    for k in range(Y.size):
        if Y[k]>y:
            j=k-1
            break
        else:
            j=k
    
    #a la interpolación sólo contribuyen los cuatro puntos de grilla más cercanos al punto (por el tamaño definido de la celda),
    #pero si está cerca del borde superior y/o derecho de la grilla, pueden ser dos o uno sólo, por eso se consideran varios
    #casos de retorno
    
    #si x,y no tocan el borde superior o derecho, contribuyen los 4 puntos
    if i!=X.size-1 and j!=Y.size-1:
        return Z.iloc[i,j]*(1-np.abs((x-X[i])/d))*(1-np.abs((y-Y[j])/d))+Z.iloc[i+1,j]*(1-np.abs((x-X[i+1])/d))*(1-np.abs((y-Y[j])/d))+Z.iloc[i,j+1]*(1-np.abs((x-X[i])/d))*(1-np.abs((y-Y[j+1])/d))+Z.iloc[i+1,j+1]*(1-np.abs((x-X[i+1])/d))*(1-np.abs((y-Y[j+1])/d))
    #si x,y tocan el borde derecho (y no superior), contribuyen los 2 puntos de arriba y abajo
    elif i==X.size-1 and j!=Y.size-1:
        return Z.iloc[i,j]*(1-np.abs((y-Y[j])/d))+Z.iloc[i,j+1]*(1-np.abs((y-Y[j+1])/d))
    #si x,y tocan el borde superior (y no derecho), contribuyen los 2 puntos de izquierda y derecha
    elif i!=X.size-1 and j==Y.size-1:
        return Z.iloc[i,j]*(1-np.abs((x-X[i])/d))+Z.iloc[i+1,j]*(1-np.abs((x-X[i+1])/d))
    #si x,y tocan ambos bordes, estará en la esquina y sólo contribuirá ese punto de grilla
    else:
        return Z.iloc[i,j]
    
#avance o push de una partícula de especie dada (ion o electrón) en un intervalo de tiempo dt, recibiendo las posiciones, 
#velocidades y factor de Lorentz x,y,u_x,u_y,gamma del paso anterior; se utilizan los campos interpolados para avanzar las
#partículas
def particle_pusher(L,gamma0,sigma_ji,species,dt,x,y,u_x,u_y,gamma,E_x,E_y,B_z):
    #la razón q/m varía entre electrones y iones
    if species=='ion':
        q_m=1
    elif species=='electron':
        q_m=-100
    
    #constante alpha que acompaña al término de la fuerza de Lorentz al adimensionalizar
    C=2*np.pi*gamma0*np.sqrt(sigma_ji)*q_m
    
    #término de la aceleración por la fuerza de Lorentz en el instante inicial
    a_x=C*(E_x+u_y*B_z/gamma)
    a_y=C*(E_y-u_x*B_z/gamma)

    #velocidades avanzadas en medio paso de tiempo usando leapfrog
    u_x_half=u_x+a_x*dt/2
    u_y_half=u_y+a_y*dt/2

    tau_z=C*dt*B_z/2
        
    u_x_prime=u_x_half+C*dt*E_x/2
    u_y_prime=u_y_half+C*dt*E_y/2
                
    sigma=1+u_x_prime**2+u_y_prime**2-tau_z**2
    
    #avance de gamma    
    gamma=np.sqrt((sigma+np.sqrt((sigma**2)+4*tau_z**2))/2)
    
    t_z=tau_z/gamma
        
    #avance de u_x,u_y en un paso entero de tiempo
    u_x=(u_x_prime+u_y_prime*t_z)/(1+t_z**2)
    u_y=(u_y_prime-u_x_prime*t_z)/(1+t_z**2)
    
    #avance de x,y en un paso entero de tiempo usando el nuevo u y gamma
    x=x+dt*u_x/gamma
    y=y+dt*u_y/gamma
    
    #condición de contorno periódica en y para partículas (L es el tamaño lineal de la región)
    if y>L/2:
        y=y-L
    elif y<-L/2:
        y=y+L
    
    #retorno de las nuevas posiciones, velocidades y factor de Lorentz en el siguiente instante de tiempo
    return x,y,u_x,u_y,gamma

#función de forma como una spline de orden 2, en una celda de tamaño 2d
def s2(x,x_p,d):
    if np.abs(x-x_p)<=d:
        return 1-np.abs((x-x_p)/d)
    else:
        return 0

#definición de la componente 1 del vector W en el algoritmo de deposición de corrientes de Esirkepov, para el punto de grilla
#X,Y, la posición antigua de la partícula x0,y0 y la posición nueva x1,y1; se usa la función de forma s2
def W1(x0,x1,y0,y1,X,Y,d):
    return 0.5*s2(x1,X,d)*s2(y1,Y,d)-0.5*s2(x0,X,d)*s2(y1,Y,d)+0.5*s2(x1,X,d)*s2(y0,Y,d)-0.5*s2(x0,X,d)*s2(y0,Y,d)

#definición de la componente 2 del vector W en el algoritmo de deposición de corrientes de Esirkepov, análoga a W1
def W2(x0,x1,y0,y1,X,Y,d):
    return 0.5*s2(x1,X,d)*s2(y1,Y,d)-0.5*s2(x1,X,d)*s2(y0,Y,d)+0.5*s2(x0,X,d)*s2(y1,Y,d)-0.5*s2(x0,X,d)*s2(y0,Y,d)

#cálculo de la suma sobre todas las partículas de una especie dada de q_p*W_pij (para cualquier componente de W) en un punto 
#de grilla fijo X,Y; recibe todas las posiciones antiguas y nuevas de las partículas
def W_total(species,component,x0,x1,y0,y1,X,Y,d):
    #la carga depende de si son iones o electrones
    if species=='ion':
        q=1
    elif species=='electron':
        q=-1
    
    #se hace la suma sobre todas las partículas dependiendo de si se quiere calcular la componente 1 (W1) o 2 (W2)
    W=0
    if component==1:
        for k in range(x0.size):
            W+=W1(x0.iloc[k],y0.iloc[k],x1.iloc[k],y1.iloc[k],X,Y,d)
        return q*W
    elif component==2:
        for k in range(x0.size):
            W+=W2(x0.iloc[k],y0.iloc[k],x1.iloc[k],y1.iloc[k],X,Y,d)
        return q*W

#algoritmo FDTD 2D implementado para avanzar los campos E_x,E_y,B_z sobre todos los puntos de grilla al siguiente instante
#de tiempo; recibe las corrientes calculadas por deposición de corriente, el paso de tiempo y espacio de grilla
def FDTD(gamma0,sigma_ji,delta,xdim,ydim,deltat,time_tot,Ex,Ey,Bz,Jx,Jy):
    
    #constante que acompaña al término de las corrientes tras adimensionalizar
    C=2*np.pi/(gamma0*np.sqrt(sigma_ji))
    
    #por simplicidad, se aplican condiciones periódicas en x y y
    
    #actualización del campo E_x en cada punto de grilla en el siguiente instante de tiempo por diferencias finitas con el
    #B_z del instante anterior, junto a la aplicación de las condiciones de contorno
    for i in range(xdim):
        for j in range(ydim-1):
            if i<xdim-1:
                Ex.iloc[i,j]=Ex.iloc[i,j]+(deltat/delta)*(Bz.iloc[i+1,j]-Bz.iloc[i,j]-Jx.iloc[i,j]*C*delta)
            elif i==xdim-1:
                Ex.iloc[i,j]=Ex.iloc[i,j]+(deltat/delta)*(Bz.iloc[i,j]-Bz.iloc[i-1,j]-Jx.iloc[i,j]*C*delta)
        Ex.iloc[i,xdim-1]=Ex.iloc[i,0]
    
    #actualización del campo E_y en cada punto de grilla en el siguiente instante de tiempo por diferencias finitas con el
    #B_z del instante anterior, junto a la aplicación de las condiciones de contorno
    for i in range(xdim):
        for j in range(ydim-1):
                Ey.iloc[i,j]=Ey.iloc[i,j]-(deltat/delta)*(Bz.iloc[i,j+1]-Bz.iloc[i,j]+Jy.iloc[i,j]*C*delta)
        Ey.iloc[i,ydim-1]=Ey.iloc[i,0]
    
    #actualización del campo B_z en cada punto de grilla en el siguiente instante de tiempo por diferencias finitas con los
    #E_x y E_y del instante nuevo, junto a la aplicación de las condiciones de contorno
    for i in range(xdim):
        for j in range(ydim-1):
            if 0<i<xdim-1 and j==0:
                Bz.iloc[i,0]=Bz.iloc[i,0]+(deltat/delta)*(Ex.iloc[i,0]-Ex.iloc[i-1,0]-Ey.iloc[i,0]+Ey.iloc[i,ydim-2])
                Bz.iloc[i,ydim-1]=Bz.iloc[i,0]
            elif i==0 and 0<j<ydim-1:
                Bz.iloc[0,j]=Bz.iloc[0,j]+(deltat/delta)*(Ex.iloc[0,j]-Ex.iloc[xdim-2,j]-Ey.iloc[0,j]+Ey.iloc[0,j-1])
                Bz.iloc[xdim-1,j]=Bz.iloc[0,j]
            elif 0<i<xdim-1 and 0<j<ydim-1:
                Bz.iloc[i,j]=Bz.iloc[i,j]+(deltat/delta)*(Ex.iloc[i,j]-Ex.iloc[i-1,j]-Ey.iloc[i,j]+Ey.iloc[i,j-1])
    
    #en las esquinas de la grilla simplemente se promedia B_z con los valores vecinos
    Bz.iloc[0,0]=(Bz.iloc[1,0]+Bz.iloc[0,1])/2
    Bz.iloc[0,ydim-1]=(Bz.iloc[0,ydim-2]+Bz.iloc[1,ydim-1])/2
    Bz.iloc[xdim-1,0]=(Bz.iloc[xdim-1,1]+Bz.iloc[xdim-2,0])/2
    Bz.iloc[xdim-1,ydim-1]=(Bz.iloc[xdim-2,ydim-1]+Bz.iloc[xdim-1,ydim-2])/2
    
    #se retornan los campos en el nuevo instante de tiempo en todo punto de grilla
    return Ex,Ey,Bz