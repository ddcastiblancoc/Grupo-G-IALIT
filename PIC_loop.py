# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 08:06:38 2021

@author: Estandar
"""

import numpy as np
import pandas as pd

import time

#importación del muestreo de u iniciales por la distribución de Maxwell-Juttner por el algoritmo de Sobol
import PIC_algorithms as PIC
#importación de los algoritmos de muestreo de la distribución de Maxwell-Juttner para las velocidades del jet,
#en específico el algoritmo de Sobol
import maxwell_juttner_algorithms as mj

#número de puntos de grilla por lado (NxN)
N=100
#largo de la grilla por lado
L=10
#espaciamiento de la grilla
dx=L/(N-1)
#intervalo de tiempo que se toma como aquel de la condición de Courant
dt=dx/np.sqrt(2)
#número de pasos de tiempo tomados
Nt=100
#anchura del jet introducido al medio
W_n=2
#factor de Lorentz de inyección del jet
gamma0=16
#velocidad de inyección del jet
beta=(1-gamma0**-2)**0.5
#temperatura adimensionalizada de los iones y electrones del jet introducido
theta_i=1
theta_e=1
#número de iones y electrones del jet introducidos al medio cada paso de tiempo
n_i_jet_dt=3
n_e_jet_dt=3
#rata de magnetización de los iones del jet
sigma_ji=1

#N puntos de grilla en X y Y cada uno (0<=x<=L, -L/2<=y<=L/2)
X = np.linspace(0,L,N)
Y = np.linspace(-L/2,L/2,N)

#importación de los archivos .csv posiciones, velocidades y campos iniciales calculados con el archivo de inicialización
particles_i_IGM = pd.read_csv('particles_i_IGM.csv')
particles_e_IGM = pd.read_csv('particles_e_IGM.csv')
E_x = pd.read_csv('E_x.csv')
E_y = pd.read_csv('E_y.csv')
B_z = pd.read_csv('B_z.csv')

columns_i_jet = ['x_i_jet','y_i_jet','u_x_i_jet','u_y_i_jet','gamma_i_jet','E_x_i_jet','E_y_i_jet','B_z_i_jet']
columns_e_jet = ['x_e_jet','y_e_jet','u_x_e_jet','u_y_e_jet','gamma_e_jet','E_x_e_jet','E_y_e_jet','B_z_e_jet']
#inicialización de DataFrames vacíos para posiciones, velocidades, factor de Lorentz y campos interpolados en las posiciones
#para las partículas del jet; estos DataFrames se llenarán con un número fijo de partículas cada paso de tiempo
particles_i_jet = pd.DataFrame([], columns=columns_i_jet)
particles_e_jet = pd.DataFrame([], columns=columns_e_jet)

#se hace una copia del campo inicial, sobre el cual se centrarán los campos en los instantes de tiempo enteros (estos se
#avanzan en pasos semienteros)
B_z_center = B_z.copy()

t_i=time.time()

#ciclo PIC que recorre Nt pasos de tiempo
for t in range(Nt):
    
    #se imprime en consola el paso de tiempo actual cuando está compilando
    print(t)
    
    #se actualizan a cero los DataFrames de densidad de corriente J_x,J_y en cada instante de tiempo
    J_x = pd.DataFrame(np.zeros([N,N]))
    J_y = pd.DataFrame(np.zeros([N,N]))
    
    #introducción de los nuevos iones y electrones del jet que serán unidos a la lista previa de partículas del jet
    particles_i_jet_new = pd.DataFrame(np.zeros([n_i_jet_dt,8]), columns=columns_i_jet)
    particles_e_jet_new = pd.DataFrame(np.zeros([n_e_jet_dt,8]), columns=columns_e_jet)
    #se inicializan las posiciones en y de las partículas del jet uniformemente en un ancho W_n alrededor de cero (se mantienen las posiciones x en cero)
    particles_i_jet_new['y_i_jet'] = np.random.uniform(-W_n/2,W_n/2,n_i_jet_dt)
    particles_e_jet_new['y_e_jet'] = np.random.uniform(-W_n/2,W_n/2,n_e_jet_dt)
    #se inicializan las velocidades en x de las partículas del jet por la distribución de Maxwell-Juttner desplazada con velocidad de inyección beta y 
    #temperatura theta (se mantienen las velocidades en y como cero; no se alteran las direcciones transversales por el boost)
    particles_i_jet_new['u_x_i_jet'] = mj.flipping_maxwell_juttner_shifted(PIC.sobol_filter(theta_i,n_i_jet_dt),beta,n_i_jet_dt)[0]
    particles_e_jet_new['u_x_e_jet'] = mj.flipping_maxwell_juttner_shifted(PIC.sobol_filter(theta_e,n_e_jet_dt),beta,n_e_jet_dt)[0]
    particles_i_jet_new['gamma_i_jet'] = np.sqrt(1+particles_i_jet_new['u_x_i_jet']**2)
    particles_e_jet_new['gamma_e_jet'] = np.sqrt(1+particles_e_jet_new['u_x_e_jet']**2)
    
    #se unen las partículas nuevas al DataFrame de partículas del jet del anterior paso
    particles_i_jet = pd.concat([particles_i_jet,particles_i_jet_new], ignore_index=True)
    particles_e_jet = pd.concat([particles_e_jet,particles_e_jet_new], ignore_index=True)
    
    #se hace una copia de las posiciones de las partículas del paso anterior (IGM y jet) para utilizarlas en la deposición de corriente
    x_i_IGM_old = particles_i_IGM['x_i_IGM'].copy()
    y_i_IGM_old = particles_i_IGM['y_i_IGM'].copy()
    x_e_IGM_old = particles_e_IGM['x_e_IGM'].copy()
    y_e_IGM_old = particles_e_IGM['y_e_IGM'].copy()
    x_i_jet_old = particles_i_jet['x_i_jet'].copy()
    y_i_jet_old = particles_i_jet['y_i_jet'].copy()
    x_e_jet_old = particles_e_jet['x_e_jet'].copy()
    y_e_jet_old = particles_e_jet['y_e_jet'].copy()
    
    #interpolación de los campos (evaluados en los puntos de grilla) del paso anterior sobre las posiciones de las partículas para usarlas
    #en el avance de las partículas; se llenan las columnas de los respectivos DataFrames; se usa el campo magnético centrado en el tiempo
    for i in particles_i_IGM.index:
        particles_i_IGM.loc[i,'E_x_i_IGM'] = PIC.interpolation(particles_i_IGM.loc[i,'x_i_IGM'],particles_i_IGM.loc[i,'y_i_IGM'],X,Y,E_x,dx)
        particles_i_IGM.loc[i,'E_y_i_IGM'] = PIC.interpolation(particles_i_IGM.loc[i,'x_i_IGM'],particles_i_IGM.loc[i,'y_i_IGM'],X,Y,E_y,dx)
        particles_i_IGM.loc[i,'B_z_i_IGM'] = PIC.interpolation(particles_i_IGM.loc[i,'x_i_IGM'],particles_i_IGM.loc[i,'y_i_IGM'],X,Y,B_z_center,dx)
        
    for i in particles_e_IGM.index:
        particles_e_IGM.loc[i,'E_x_e_IGM'] = PIC.interpolation(particles_e_IGM.loc[i,'x_e_IGM'],particles_e_IGM.loc[i,'y_e_IGM'],X,Y,E_x,dx)
        particles_e_IGM.loc[i,'E_y_e_IGM'] = PIC.interpolation(particles_e_IGM.loc[i,'x_e_IGM'],particles_e_IGM.loc[i,'y_e_IGM'],X,Y,E_y,dx)
        particles_e_IGM.loc[i,'B_z_e_IGM'] = PIC.interpolation(particles_e_IGM.loc[i,'x_e_IGM'],particles_e_IGM.loc[i,'y_e_IGM'],X,Y,B_z_center,dx)
    
    for i in particles_i_jet.index:
        particles_i_jet.loc[i,'E_x_i_jet'] = PIC.interpolation(particles_i_jet.loc[i,'x_i_jet'],particles_i_jet.loc[i,'y_i_jet'],X,Y,E_x,dx)
        particles_i_jet.loc[i,'E_y_i_jet'] = PIC.interpolation(particles_i_jet.loc[i,'x_i_jet'],particles_i_jet.loc[i,'y_i_jet'],X,Y,E_y,dx)
        particles_i_jet.loc[i,'B_z_i_jet'] = PIC.interpolation(particles_i_jet.loc[i,'x_i_jet'],particles_i_jet.loc[i,'y_i_jet'],X,Y,B_z_center,dx)
        
    for i in particles_e_jet.index:
        particles_e_jet.loc[i,'E_x_e_jet'] = PIC.interpolation(particles_e_jet.loc[i,'x_e_jet'],particles_e_jet.loc[i,'y_e_jet'],X,Y,E_x,dx)
        particles_e_jet.loc[i,'E_y_e_jet'] = PIC.interpolation(particles_e_jet.loc[i,'x_e_jet'],particles_e_jet.loc[i,'y_e_jet'],X,Y,E_y,dx)
        particles_e_jet.loc[i,'B_z_e_jet'] = PIC.interpolation(particles_e_jet.loc[i,'x_e_jet'],particles_e_jet.loc[i,'y_e_jet'],X,Y,B_z_center,dx)
    
    #avance de las partículas actualizando las nuevas posiciones, velocidades y factores de Lorentz con los campos interpolados
    for i in particles_i_IGM.index:
        particles_i_IGM.loc[i,'x_i_IGM'],particles_i_IGM.loc[i,'y_i_IGM'],particles_i_IGM.loc[i,'u_x_i_IGM'],particles_i_IGM.loc[i,'u_y_i_IGM'],particles_i_IGM.loc[i,'gamma_i_IGM'] = PIC.particle_pusher(L,gamma0,sigma_ji,'ion',dt,particles_i_IGM.loc[i,'x_i_IGM'],particles_i_IGM.loc[i,'y_i_IGM'],particles_i_IGM.loc[i,'u_x_i_IGM'],particles_i_IGM.loc[i,'u_y_i_IGM'],particles_i_IGM.loc[i,'gamma_i_IGM'],particles_i_IGM.loc[i,'E_x_i_IGM'],particles_i_IGM.loc[i,'E_y_i_IGM'],particles_i_IGM.loc[i,'B_z_i_IGM'])
        
    for i in particles_e_IGM.index:
        particles_e_IGM.loc[i,'x_e_IGM'],particles_e_IGM.loc[i,'y_e_IGM'],particles_e_IGM.loc[i,'u_x_e_IGM'],particles_e_IGM.loc[i,'u_y_e_IGM'],particles_e_IGM.loc[i,'gamma_e_IGM'] = PIC.particle_pusher(L,gamma0,sigma_ji,'electron',dt,particles_e_IGM.loc[i,'x_e_IGM'],particles_e_IGM.loc[i,'y_e_IGM'],particles_e_IGM.loc[i,'u_x_e_IGM'],particles_e_IGM.loc[i,'u_y_e_IGM'],particles_e_IGM.loc[i,'gamma_e_IGM'],particles_e_IGM.loc[i,'E_x_e_IGM'],particles_e_IGM.loc[i,'E_y_e_IGM'],particles_e_IGM.loc[i,'B_z_e_IGM'])
    
    for i in particles_i_jet.index:
        particles_i_jet.loc[i,'x_i_jet'],particles_i_jet.loc[i,'y_i_jet'],particles_i_jet.loc[i,'u_x_i_jet'],particles_i_jet.loc[i,'u_y_i_jet'],particles_i_jet.loc[i,'gamma_i_jet'] = PIC.particle_pusher(L,gamma0,sigma_ji,'ion',dt,particles_i_jet.loc[i,'x_i_jet'],particles_i_jet.loc[i,'y_i_jet'],particles_i_jet.loc[i,'u_x_i_jet'],particles_i_jet.loc[i,'u_y_i_jet'],particles_i_jet.loc[i,'gamma_i_jet'],particles_i_jet.loc[i,'E_x_i_jet'],particles_i_jet.loc[i,'E_y_i_jet'],particles_i_jet.loc[i,'B_z_i_jet'])
        
    for i in particles_e_jet.index:
        particles_e_jet.loc[i,'x_e_jet'],particles_e_jet.loc[i,'y_e_jet'],particles_e_jet.loc[i,'u_x_e_jet'],particles_e_jet.loc[i,'u_y_e_jet'],particles_e_jet.loc[i,'gamma_e_jet'] = PIC.particle_pusher(L,gamma0,sigma_ji,'electron',dt,particles_e_jet.loc[i,'x_e_jet'],particles_e_jet.loc[i,'y_e_jet'],particles_e_jet.loc[i,'u_x_e_jet'],particles_e_jet.loc[i,'u_y_e_jet'],particles_e_jet.loc[i,'gamma_e_jet'],particles_e_jet.loc[i,'E_x_e_jet'],particles_e_jet.loc[i,'E_y_e_jet'],particles_e_jet.loc[i,'B_z_e_jet'])

    #condición de contorno de partículas para x: se escogen las partículas que después del avance siguen dentro de la región, las que se salieron son descartadas
    particles_i_IGM = particles_i_IGM[np.logical_and(particles_i_IGM['x_i_IGM']>=0,particles_i_IGM['x_i_IGM']<=L)]
    particles_e_IGM = particles_e_IGM[np.logical_and(particles_e_IGM['x_e_IGM']>=0,particles_e_IGM['x_e_IGM']<=L)]
    particles_i_jet = particles_i_jet[np.logical_and(particles_i_jet['x_i_jet']>=0,particles_i_jet['x_i_jet']<=L)]
    particles_e_jet = particles_e_jet[np.logical_and(particles_e_jet['x_e_jet']>=0,particles_e_jet['x_e_jet']<=L)]
    
    #después de filtrar las partículas que siguen dentro de la región, se actualizan las posiciones antiguas asociadas a las nuevas que pasaron el filtro 
    #para que las listas de partículas tengan el mismo tamaño
    x_i_IGM_old = x_i_IGM_old[x_i_IGM_old.index.isin(particles_i_IGM.index)]
    y_i_IGM_old = y_i_IGM_old[y_i_IGM_old.index.isin(particles_i_IGM.index)]
    x_e_IGM_old = x_e_IGM_old[x_e_IGM_old.index.isin(particles_e_IGM.index)]
    y_e_IGM_old = y_e_IGM_old[y_e_IGM_old.index.isin(particles_e_IGM.index)]
    x_i_jet_old = x_i_jet_old[x_i_jet_old.index.isin(particles_i_jet.index)]
    y_i_jet_old = y_i_jet_old[y_i_jet_old.index.isin(particles_i_jet.index)]
    x_e_jet_old = x_e_jet_old[x_e_jet_old.index.isin(particles_e_jet.index)]
    y_e_jet_old = y_e_jet_old[y_e_jet_old.index.isin(particles_e_jet.index)]
    
    #se llenan los DataFrames de densidad de corriente en cada punto de grilla por las 4 contribuciones de i-IGM, e-IGM, i-jet, e-jet, por medio de las componentes de W
    for i in range(N-1):
        for j in range(N):
            J_x.iloc[i+1,j]=J_x.iloc[i,j]-(dx/dt)*(PIC.W_total('ion',1,x_i_IGM_old,particles_i_IGM['x_i_IGM'],y_i_IGM_old,particles_i_IGM['y_i_IGM'],X[i],Y[j],dx)+PIC.W_total('electron',1,x_e_IGM_old,particles_e_IGM['x_e_IGM'],y_e_IGM_old,particles_e_IGM['y_e_IGM'],X[i],Y[j],dx)+PIC.W_total('ion',1,x_i_jet_old,particles_i_jet['x_i_jet'],y_i_jet_old,particles_i_jet['y_i_jet'],X[i],Y[j],dx)+PIC.W_total('electron',1,x_e_jet_old,particles_e_jet['x_e_jet'],y_e_jet_old,particles_e_jet['y_e_jet'],X[i],Y[j],dx))

    for i in range(N):
        for j in range(N-1):
            J_y.iloc[i,j+1]=J_x.iloc[i,j]-(dx/dt)*(PIC.W_total('ion',2,x_i_IGM_old,particles_i_IGM['x_i_IGM'],y_i_IGM_old,particles_i_IGM['y_i_IGM'],X[i],Y[j],dx)+PIC.W_total('electron',2,x_e_IGM_old,particles_e_IGM['x_e_IGM'],y_e_IGM_old,particles_e_IGM['y_e_IGM'],X[i],Y[j],dx)+PIC.W_total('ion',2,x_i_jet_old,particles_i_jet['x_i_jet'],y_i_jet_old,particles_i_jet['y_i_jet'],X[i],Y[j],dx)+PIC.W_total('electron',2,x_e_jet_old,particles_e_jet['x_e_jet'],y_e_jet_old,particles_e_jet['y_e_jet'],X[i],Y[j],dx))
    
    #antes de actualizar el campo magnético al siguiente instante de tiempo, se guarda una copia para poder centrar el campo en el paso entero siguiente
    B_z_old = B_z.copy()
    
    #se avanzan los campos eléctrico y magnético por FDTD recibiendo los campos del paso anterior y las corrientes de deposición de corriente
    E_x,E_y,B_z = PIC.FDTD(gamma0,sigma_ji,dx,N,N,dt,Nt,E_x,E_y,B_z,J_x,J_y)
    
    #se centra el campo magnético en el siguiente paso entero; se guardan el campo magnético de pasos semienteros B_z y de pasos enteros B_z_center
    #para actualizarlo con FDTD y usar el pusher respectivamente
    B_z_center = (B_z_old+B_z)/2
    
print(time.time()-t_i) 

#una vez acaba el ciclo PIC en un tiempo Nt*dt, se exportan los DataFrames resultantes de partículas y campos para cálculos/gráficas posteriores
particles_i_IGM.to_csv('particles_i_IGM_final.csv', index=False)
particles_e_IGM.to_csv('particles_e_IGM_final.csv', index=False)
particles_i_jet.to_csv('particles_i_jet_final.csv', index=False)
particles_e_jet.to_csv('particles_e_jet_final.csv', index=False)
E_x.to_csv('E_x_final.csv', index=False)
E_y.to_csv('E_y_final.csv', index=False)
B_z.to_csv('B_z_final.csv', index=False)