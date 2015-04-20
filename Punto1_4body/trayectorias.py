import numpy as np
import matplotlib.pyplot as plt

print "Importando datos"
#Importar datos
datos = np.genfromtxt('orbitas.txt')

print "Creando grafica"
plt.plot(datos[:,0], datos[:,1], color = 'r', label= 'Sol')
plt.plot(datos[:,3], datos[:,4], color = 'g', label = 'Tierra')
plt.plot(datos[:,6], datos[:,7], color = 'b', label = 'Luna')
plt.plot(datos[:,9], datos[:,10], color = 'y', label = 'Asteroide')

plt.legend()

plt.show()


