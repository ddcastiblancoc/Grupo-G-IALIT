# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 14:31:47 2021

@author: Estandar
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#texto TeX en gráficas (des-comentar si tiene un compilador TeX instalado, p. ej., Texmaker)
"""
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": "Computer Modern Roman"})
"""

#número de puntos de grilla por lado (NxN)
N=100
#largo de la grilla por lado
L=10
#N puntos de grilla en X y Y cada uno (0<=x<=L, -L/2<=y<=L/2)
X = np.linspace(0,L,N)
Y = np.linspace(-L/2,L/2,N)

#-------------------------------------------------------------------------------------------------------------------------------------------------------------
#si compiló los códigos anteriores de cero y obtuvo los archivos .csv finales con PIC_loop.py, esta sección dará las respectivas gráficas
#-------------------------------------------------------------------------------------------------------------------------------------------------------------

#importación de los archivos finales .csv del ciclo PIC de posiciones, velocidades y campos iniciales calculados con el archivo 
#de inicialización
particles_i_IGM_final = pd.read_csv('particles_i_IGM_final.csv')
particles_e_IGM_final = pd.read_csv('particles_e_IGM_final.csv')
particles_i_jet_final = pd.read_csv('particles_i_jet_final.csv')
particles_e_jet_final = pd.read_csv('particles_e_jet_final.csv')
E_x_final = pd.read_csv('E_x_final.csv')
E_y_final = pd.read_csv('E_y_final.csv')
B_z_final = pd.read_csv('B_z_final.csv')

#gráfica de (1) las posiciones finales del sistema del jet (rojo) y IGM (azul), y el contorno del campo magnético final
plt.figure()

plt.subplot(2,1,1)
plt.scatter(particles_i_IGM_final['x_i_IGM'],particles_i_IGM_final['y_i_IGM'], c=particles_i_IGM_final['B_z_i_IGM'], cmap='Blues')
plt.scatter(particles_e_IGM_final['x_e_IGM'],particles_e_IGM_final['y_e_IGM'], c=particles_e_IGM_final['B_z_e_IGM'], cmap='Blues')
plt.scatter(particles_i_jet_final['x_i_jet'],particles_i_jet_final['y_i_jet'], c=particles_i_jet_final['B_z_i_jet'], cmap='Reds')
plt.scatter(particles_e_jet_final['x_e_jet'],particles_e_jet_final['y_e_jet'], c=particles_e_jet_final['B_z_e_jet'], cmap='Reds')
plt.xlabel(r'$x/\lambda_{ji}$')
plt.ylabel(r'$y/\lambda_{ji}$')

plt.subplot(2,1,2)
plt.contourf(X,Y,B_z_final.T, 20,cmap='plasma')
plt.colorbar()
plt.xlabel(r'$x/\lambda_{ji}$')
plt.ylabel(r'$y/\lambda_{ji}$')

#-------------------------------------------------------------------------------------------------------------------------------------------------------------
#si quiere compilar los códigos con los datos que nosotros obtuvimos de la carpeta dm, esta sección dará las respectivas gráficas si des-comenta esta sección
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
#importación de los archivos finales .csv del ciclo PIC de posiciones, velocidades y campos iniciales calculados con el archivo 
#de inicialización
particles_i_IGM_final = pd.read_csv('dm/particles_i_IGM_final.csv')
particles_e_IGM_final = pd.read_csv('dm/particles_e_IGM_final.csv')
particles_i_jet_final = pd.read_csv('dm/particles_i_jet_final.csv')
particles_e_jet_final = pd.read_csv('dm/particles_e_jet_final.csv')
E_x_final = pd.read_csv('dm/E_x_final.csv')
E_y_final = pd.read_csv('dm/E_y_final.csv')
B_z_final = pd.read_csv('dm/B_z_final.csv')

#gráfica de (1) las posiciones finales del sistema del jet (rojo) y IGM (azul), y el contorno del campo magnético final
plt.figure()

plt.subplot(2,1,1)
plt.scatter(particles_i_IGM_final['x_i_IGM'],particles_i_IGM_final['y_i_IGM'], c=particles_i_IGM_final['B_z_i_IGM'], cmap='Blues')
plt.scatter(particles_e_IGM_final['x_e_IGM'],particles_e_IGM_final['y_e_IGM'], c=particles_e_IGM_final['B_z_e_IGM'], cmap='Blues')
plt.scatter(particles_i_jet_final['x_i_jet'],particles_i_jet_final['y_i_jet'], c=particles_i_jet_final['B_z_i_jet'], cmap='Reds')
plt.scatter(particles_e_jet_final['x_e_jet'],particles_e_jet_final['y_e_jet'], c=particles_e_jet_final['B_z_e_jet'], cmap='Reds')
plt.xlabel(r'$x/\lambda_{ji}$')
plt.ylabel(r'$y/\lambda_{ji}$')

plt.subplot(2,1,2)
plt.contourf(X,Y,B_z_final.T, 20,cmap='plasma')
plt.colorbar()
plt.xlabel(r'$x/\lambda_{ji}$')
plt.ylabel(r'$y/\lambda_{ji}$')
"""
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
#si quiere compilar los códigos con los datos que nosotros obtuvimos de la carpeta fm, esta sección dará las respectivas gráficas si des-comenta esta sección
#-------------------------------------------------------------------------------------------------------------------------------------------------------------
"""
#importación de los archivos finales .csv del ciclo PIC de posiciones, velocidades y campos iniciales calculados con el archivo 
#de inicialización
particles_i_IGM_final = pd.read_csv('fm/particles_i_IGM_final.csv')
particles_e_IGM_final = pd.read_csv('fm/particles_e_IGM_final.csv')
particles_i_jet_final = pd.read_csv('fm/particles_i_jet_final.csv')
particles_e_jet_final = pd.read_csv('fm/particles_e_jet_final.csv')
E_x_final = pd.read_csv('fm/E_x_final.csv')
E_y_final = pd.read_csv('fm/E_y_final.csv')
B_z_final = pd.read_csv('fm/B_z_final.csv')

#gráfica de (1) las posiciones finales del sistema del jet (rojo) y IGM (azul), y el contorno del campo magnético final
plt.figure()

plt.subplot(2,1,1)
plt.scatter(particles_i_IGM_final['x_i_IGM'],particles_i_IGM_final['y_i_IGM'], c=particles_i_IGM_final['B_z_i_IGM'], cmap='Blues')
plt.scatter(particles_e_IGM_final['x_e_IGM'],particles_e_IGM_final['y_e_IGM'], c=particles_e_IGM_final['B_z_e_IGM'], cmap='Blues')
plt.scatter(particles_i_jet_final['x_i_jet'],particles_i_jet_final['y_i_jet'], c=particles_i_jet_final['B_z_i_jet'], cmap='Reds')
plt.scatter(particles_e_jet_final['x_e_jet'],particles_e_jet_final['y_e_jet'], c=particles_e_jet_final['B_z_e_jet'], cmap='Reds')
plt.xlabel(r'$x/\lambda_{ji}$')
plt.ylabel(r'$y/\lambda_{ji}$')

plt.subplot(2,1,2)
plt.contourf(X,Y,B_z_final.T, 20,cmap='plasma')
plt.colorbar()
plt.xlabel(r'$x/\lambda_{ji}$')
plt.ylabel(r'$y/\lambda_{ji}$')
"""