# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 22:54:50 2021

@author: Estandar
"""

import numpy as np
import pandas as pd

#importación de las funciones y algoritmos a utilizar en la inicialización y el ciclo PIC
import PIC_algorithms as PIC

#número de puntos de grilla por lado (NxN)
N=100
#largo de la grilla por lado
L=10
#espaciamiento de la grilla
dx=L/(N-1)
#número de iones iniciales en el medio intergaláctico magnetizado (IGM)
n_initial_i_IGM=50
#número de electrones iniciales en el IGM
n_initial_e_IGM=50
#ancho de la región de cambio de polaridad del campo magnético inicial
W_b=1

#N puntos de grilla en X y Y cada uno (0<=x<=L, -L/2<=y<=L/2)
X = np.linspace(0,L,N)
Y = np.linspace(-L/2,L/2,N)

columns_i_IGM = ['x_i_IGM','y_i_IGM','u_x_i_IGM','u_y_i_IGM','gamma_i_IGM','E_x_i_IGM','E_y_i_IGM','B_z_i_IGM']
columns_e_IGM = ['x_e_IGM','y_e_IGM','u_x_e_IGM','u_y_e_IGM','gamma_e_IGM','E_x_e_IGM','E_y_e_IGM','B_z_e_IGM']

#DataFrames de Pandas que almacenan las posiciones, velocidades, factor de Lorentz iniciales de iones y electrones 
#del medio (por separado); se asignan columnas adicionales para posteriormente llenar con los campos E_x,E_y,B_z interpolados
#en la posición de las partículas; estos DataFrames se actualizan en cada instante de tiempo del ciclo PIC
particles_i_IGM = pd.DataFrame(np.zeros([n_initial_i_IGM,8]), columns=columns_i_IGM)
particles_e_IGM = pd.DataFrame(np.zeros([n_initial_e_IGM,8]), columns=columns_e_IGM)

#se asignan DataFrames con columnas genéricas para cada punto de grilla para los campos
E_x = pd.DataFrame(np.zeros([N,N]))
E_y = pd.DataFrame(np.zeros([N,N]))
B_z = pd.DataFrame(np.zeros([N,N]))

#se asigma un DataFrame para la densidad de carga inicial con el que hallar el campo eléctrico inicial
rho = pd.DataFrame(np.zeros([N,N]))

#se inicializan las posiciones de las partículas en el medio uniformemente; las velocidades se mantienen en cero (por ende 
#los factores de Lorentz son todos 1)
particles_i_IGM['x_i_IGM'] = np.random.uniform(0,L,n_initial_i_IGM)
particles_i_IGM['y_i_IGM'] = np.random.uniform(-L/2,L/2,n_initial_i_IGM)
particles_e_IGM['x_e_IGM'] = np.random.uniform(0,L,n_initial_e_IGM)
particles_e_IGM['y_e_IGM'] = np.random.uniform(-L/2,L/2,n_initial_e_IGM)
particles_i_IGM['gamma_i_IGM'] = np.ones(n_initial_i_IGM)
particles_e_IGM['gamma_e_IGM'] = np.ones(n_initial_e_IGM)

#se juntan las posiciones x de todas las partículas en una sola columna, igual con y, para usarlas más fácilmente en la
#proyección de densidad de carga en la grilla
x_IGM = pd.concat([particles_i_IGM['x_i_IGM'],particles_e_IGM['x_e_IGM']], ignore_index=True)
y_IGM = pd.concat([particles_i_IGM['y_i_IGM'],particles_e_IGM['y_e_IGM']], ignore_index=True)

#aplicación de la proyección de densidad de carga en la grilla utilizando las posiciones de las partículas del medio definidas
#anteriormente, sobre cada punto de la grilla X,Y; el tamaño de las celdas se toma como el paso de la grilla
for i in range(N):
    for j in range(N):
        rho.iloc[i,j] = PIC.projection(X[i],Y[j],x_IGM,y_IGM,dx)

#se calcula por sobre-relajación el campo eléctrico inicial con la densidad de carga antes calculada
E_x,E_y = PIC.overrelaxation(N,dx,X,Y,rho)
#se inicializa el campo magnético con el perfil de alternación de polaridad requerido 
B_z = PIC.B_z_initial(X,Y,W_b)

#se extraen los DataFrames de partículas y campos a archivos externos .csv para utilizarlos en el ciclo PIC por aparte
particles_i_IGM.to_csv('particles_i_IGM.csv', index=False)
particles_e_IGM.to_csv('particles_e_IGM.csv', index=False)
E_x.to_csv('E_x.csv', index=False)
E_y.to_csv('E_y.csv', index=False)
B_z.to_csv('B_z.csv', index=False)