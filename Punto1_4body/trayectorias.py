#!/usr/bin/env python
# -*- coding: 850-*-

import numpy as np
import matplotlib.pyplot as plt

print "Importando datos"
#Importar datos
datos1 = np.genfromtxt('orbitas_1year.txt')
datos1000 = np.genfromtxt('orbitas_1000year.txt')

print "Creando grafica para 1 year"
fig = plt.figure()
plt.plot(datos1[:,0], datos1[:,1], color = 'r', label= 'Sol')
plt.plot(datos1[:,3], datos1[:,4], color = 'g', label = 'Tierra')
plt.plot(datos1[:,6], datos1[:,7], color = 'b', label = 'Luna')
plt.plot(datos1[:,9], datos1[:,10], color = 'y', label = 'Asteroide')
plt.title("Trayectorias para 1 year terrestre")
plt.legend()
plt.xlabel('Posicion en x')
plt.ylabel('Posicion en y')
plt.savefig("orbitas-1yr.png")

print "Creando grafica para 1000 year"
fig = plt.figure()
plt.plot(datos1000[:,0], datos1000[:,1], color = 'r', label= 'Sol')
plt.plot(datos1000[:,3], datos1000[:,4], color = 'g', label = 'Tierra')
plt.plot(datos1000[:,6], datos1000[:,7], color = 'b', label = 'Luna')
plt.plot(datos1000[:,9], datos1000[:,10], color = 'y', label = 'Asteroide')
plt.title("Trayectorias para 1000 years terrestres")
plt.xlabel('Posicion en x')
plt.ylabel('Posicion en y')
plt.legend()

plt.savefig("orbitas-1000yr.png")



