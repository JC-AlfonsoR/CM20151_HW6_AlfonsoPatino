#!/usr/bin/env python
# -*- coding: 850-*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

print "Importando datos"
#Importar datos
datos1 = np.genfromtxt('orbitas_1year.txt')
datos1000 = np.genfromtxt('orbitas_1000year.txt')

print "Creando grafica para 1 year"
fig = plt.figure()
plt.scatter(datos1[:,0], datos1[:,1], color = 'r', label= 'Sol', s=1)
plt.scatter(datos1[:,3], datos1[:,4], color = 'g', label = 'Tierra', s=1)
plt.scatter(datos1[:,6], datos1[:,7], color = 'b', label = 'Luna', s=1)
plt.scatter(datos1[:,9], datos1[:,10], color = 'y', label = 'Asteroide', s=1)
plt.title("Trayectorias para 1 year terrestre")
plt.legend()
plt.xlabel('Posicion en x (Au)')
plt.ylabel('Posicion en y (Au)')
plt.savefig("orbitas-1yr.png")

print "Creando grafica para 1000 year"
fig = plt.figure()
plt.scatter(datos1000[:,0], datos1000[:,1], color = 'r', label= 'Sol', s=1)
plt.scatter(datos1000[:,3], datos1000[:,4], color = 'g', label = 'Tierra', s=1)
plt.scatter(datos1000[:,6], datos1000[:,7], color = 'b', label = 'Luna', s=1)
plt.scatter(datos1000[:,9], datos1000[:,10], color = 'y', label = 'Asteroide', s=1)
plt.title("Trayectorias para 1000 years terrestres")
plt.xlabel('Posicion en x (Au)')
plt.ylabel('Posicion en y (Au)')
plt.legend()

plt.savefig("orbitas-1000yr.png")






