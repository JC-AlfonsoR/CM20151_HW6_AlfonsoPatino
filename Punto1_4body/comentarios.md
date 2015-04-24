
La solución en C toma como primer argumento el nombre del archivo en el que se encuentran las condiciones iniciales que debe estar en la misma carpeta del
ejecutable que se genera. El segundo argumento es el paso en el tiempo que se quiere usar para hacer el cálculo. El tercer argumento
corresponde al numero de años para el que se quiere hacer el calculo. Como en el enunciado piden en ambos casos un numeros entero de años, el código
no admite numero de añoss que no sean enteros. Tanto el el paso de tiempo como el numero de años deben estar en unidades de años.

El código en C está escrito en términos de las siguientes unidades:
	masa: masa del sol
	distancia: Au
	tiempo: años


Como en la solución de C usa la función apppend para escribir el archivo de las posiciones de los cuerpos, en el makefile se incluýó una función
clean que se debe ejecutar antes de ejecutar toda la secuencia de reglas usando make. De esta forma se evita que en los arvhivos de las posiciones queden
consignadas dos veces las posiciones.

Debido a la gran velocidad inicial del asteroide, en la gráfica que se genera no se puede apreciar correctamente el movimiento de la Tierra y la Luna. Sin embargo,
si se modifican los límites del gráfico se puede apreciar que las trayectorias tanto de la Tierra y la Luna son las esperadas.


