#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//Funcion que importa los datos del archivo con nombre que se le pasa como parametro.
double * importacion_datos(char *filename);


/*
Funcion que calcula la aceleracion en cada una de la tres componentes para el cuerpo 0.
Inputs:
  m_1: masa del cuerpo 1
  m_2: masa del cuerpo 2
  m_3: masa del cuerpo 3
  p0: vector posicion del cuerpo al que se le quiere calcular la aceleracion
  p1: vector posicion del cuerpo 1
  p2: vector posicion del cuerpo 2
  p3: vector posicion del cuerpo 3
*/
double *ecuacion_gravedad(double m_1, double m_2, double m_3, double *p_0, double *p_1, double *p_2, double *p_3);

//Implementacion de Runge_Kutta4 para el un time-step pequeno en una dimension de un solo cuerpo
double *rungekutta4(double p_cuerpo, double v_cuerpo);

//Calculo de funcion de la derivada de la posicion en cada una de las componentes para llamarla en la funcion rungekutta.
double *ecuacion_velocidad(double v_x_cuerpo, double v_y_cuerpo, double v_z_cuerpo);

int main(int argc, char **argv){

  char *filename;
  double *datos;
  double x_sol, y_sol, z_sol, x_tierra, y_tierra, z_tierra, x_luna, y_luna, z_luna, x_ast, y_ast, z_ast;
  double v_x_sol, v_y_sol, v_z_sol, v_x_tierra, v_y_tierra, v_z_tierra, v_x_luna, v_y_luna, v_z_luna, v_x_ast, v_y_ast, v_z_ast;
  double masa_sol, masa_luna, masa_tierra, masa_ast;

  //Memory allocation de los punteros
  filename = malloc(10*sizeof(char));

  //Inicializacion de variables
  filename = argv[1];

  //Se llama la funcion para importar los datos. Mirar el archivo ic.txt para mirar a que corresponde cada valor en el arreglo
  datos  = importacion_datos(filename);

  //Inicializacion de condiciones iniciales
  x_sol = datos[0];
  y_sol = datos[1];
  z_sol = datos[2];
  x_tierra = datos[3];
  y_tierra = datos[4];
  z_tierra = datos[5];
  x_luna = datos[6];
  y_luna = datos[7];
  z_luna = datos[8];
  x_ast = datos[9];
  y_ast = datos[10];
  z_ast = datos[11];
  v_x_sol = datos[12];
  v_y_sol = datos[13];
  v_z_sol = datos[14];
  v_x_tierra = datos[15];
  v_y_tierra = datos[16];
  v_z_tierra = datos[17];
  v_x_luna = datos[18];
  v_y_luna = datos[19];
  v_z_luna = datos[20];
  v_x_ast = datos[21];
  v_y_ast = datos[22];
  v_z_ast = datos[23];
  masa_sol = datos[24];
  masa_tierra = datos[25];
  masa_luna = datos[26];
  masa_ast = datos[27];


  return 0;

}


double * importacion_datos(char *filename){

  FILE *in;
  int i;
  char var[50];
  double *datos;
  char valor[10]; 

  //Memory allocation de los punteros
  datos = malloc(27*sizeof(double));  

  in = fopen(filename, "r");
  
  if(!in)
    {
      printf("Problemas leyendo archivo %s \n", filename);
      exit(1);
      
    }

  //Se lee el archivo para importar los datos
  for(i=0;i<28;i++)
    {
      fscanf(in, "%s \t %s \n", &var, &valor);
      datos[i] = atof(valor);

    }

  fclose(in); 
  
  return datos;
}

double *ecuacion_gravedad(double m_1, double m_2, double m_3, double *p0, double *p1, double *p2, double *p3){

  double *dydt;
  double G, denominador1, denominador2, denominador3; //Constante de gravitacion universal


  dydt = malloc(3*sizeof(double));

  G = 6.67384E-11;

  denominador1 = sqrt(pow(pow(p0[0] - p1[0], 2 ) + pow ( p0[1] - p1[1], 2 ) + pow(p0[2] - p1[2] , 2), 3 ) );
  denominador2 = sqrt(pow(pow(p0[0] - p2[0], 2 ) + pow ( p0[1] - p2[1], 2 ) + pow(p0[2] - p2[2] , 2), 3 ) );
  denominador3 = sqrt(pow(pow(p0[0] - p3[0], 2 ) + pow ( p0[1] - p3[1], 2 ) + pow(p0[2] - p3[2] , 2), 3 ) );

  //Calculo para componente en x
  dydt[0] = G*((m_1*(p0[0]-p1[0])/denominador1) + (m_2*(p0[0]-p2[0])/denominador2) + (m_3*(p0[0]-p3[0])/denominador3));

  //Calculo para componente en y
  dydt[1] = G*((m_1*(p0[1]-p1[1])/denominador1) + (m_2*(p0[1]-p2[1])/denominador2) + (m_3*(p0[1]-p3[1])/denominador3));

  //Calculo para componente en z
  dydt[2] = G*((m_1*(p0[2]-p1[2])/denominador1) + (m_2*(p0[2]-p2[2])/denominador2) + (m_3*(p0[2]-p3[2])/denominador3));

  return dydt;
}

double *ecuacion_velocidad(double v_x_cuerpo, double v_y_cuerpo, double v_z_cuerpo){

  double *dydt;

  dydt = malloc(3*sizeof(double));

  dydt[0] = v_x_cuerpo;
  dydt[1] = v_y_cuerpo;
  dydt[2] = v_z_cuerpo;  

  return dydt;
}

/*
Se debe tener en cuenta la posicion de la particula en la componente
y la distancia (denominador) para cada una de las particulas. Se escribe
por ahora en terminos de los denominadores. Se debe modificar la funcion
ecuacion_gravedad
 */

double *rungekutta4(double p_cuerpo, double den1, double den2, double den3, double v_cuerpo, double m1, double m2, double m3, double tiempo){

  double *nuevas_coordenadas;
  double k1_posicion, k1_velocidad, k2_posicion, k2_velocidad,k3_posicion, k3_velocidad,k4_posicion, k4_velocidad;

  nuevas_coordenadas = malloc(3*sizeof(double));

  k1_posicion = ecuacion_velocidad(v_cuerpo);
  k1_velocidad = ecuacion_gravedad(m1, m2, m3, p_cuerpo, den1, den2, den3);

  //Primer paso
  tiempo_1 = tiempo + h*0.5;
  p_1 = p_cuerpo + h*0.5*k1_posicion;
  v_1 = v_cuerpo + h*0.5*k1_velocidad;
  k2_posicion = ecuacion_velocidad(v_1);
  k2_velocidad = ecuacion_gravedad(m1, m2, m3, p_1, den1, den2, den3);

  //Segundo paso
  tiempo_2 = tiempo + h*0.5;
  p_2 = p_cuerpo + h*0.5*k2_posicion;
  v_2 = v_cuerpo + h*0.5*k2_velocidad;
  k3_posicion = ecuacion_velocidad(v_2);
  k3_velocidad = ecuacion_gravedad(m1, m2, m3, p_2, den1, den2, den3);

  //Tercer paso
  tiempo_3 = tiempo + h;
  p_3 = p_cuerpo + h*k3_posicion;
  v_3 = v_cuerpo + h*k3_velocidad;
  k4_posicion = ecuacion_velocidad(v_1);
  k4_velocidad = ecuacion_gravedad(m1, m2, m3, p_3, den1, den2, den3);

  //Cuarto paso
  kprom_posicion = (1.0/6.0)*(k1_posicion+2.0*k2_posicion+2.0*k3_posicion+k4_posicion);
  kprom_velocidad = (1.0/6.0)*(k1_velocidad+2.0*k2_velocidad+2.0*k3_velocidad+k4_velocidad);

  //Creacion de nuevos valores
  p_nueva = p_cuerpo + h*kprom_posicion;
  v_nueva = v_cuerpo + v_cuerpo + h*kprom_velocidad;

  nuevas_coordenadas[0] = p_nueva;
  nuevas_coordenadas[1] = v_nueva;

  return nuevas_coordenadas;

} 




  


