#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define G_GRAV 39.486//6.67384E-11;
#define FLOAT float

/*
Se adopta la siguiente convencion de indices. 0: es para el cuerpo que se mueve. 1,2,3 son para los cuerpos que estan afectando al cuerpo
en movimiento
*/

/*
Funcion que importa los datos desde el archivo de condiciones iniciales.
Input: nombre del archivo
Output: array con condiciones iniciales. (Para ver a que corresponde cada valor en el array ver archivo ic.txt)
*/
FLOAT * importacion_datos(char *filename);

/*
Funcion que calcula la distancia entre el cuerpo 0 con cada una de las particulas. Teniendo este valor se calcula la norma de esta distancia
y se eleva a la tres. Esto corresponderia a los denominadores de la ecuacion de gravedad.
Inputs:
  p0: vector posicion del cuerpo 0
  p1: vector posicion del cuerpo 1
  p2: vector posicion del cuerpo 2
  p3: vector posicion del cuerpo 3
Output:
  Los tres denominadores de la ecuacion de gravedad.
*/
 FLOAT *calcular_distancia(FLOAT x0, FLOAT y0, FLOAT z0, FLOAT x1, FLOAT y1, FLOAT z1, FLOAT x2, FLOAT y2, FLOAT z2, FLOAT x3, FLOAT y3, FLOAT z3);

/*
Funcion que calcula la aceleracion en una componente para el cuerpo 0.
Inputs:
  m_1: masa del cuerpo 1
  m_2: masa del cuerpo 2
  m_3: masa del cuerpo 3
  p0: posicion del cuerpo 0
  p1: posicion del cuerpo 1
  p2: posicion del cuerpo 2
  p3: posicion del cuerpo 3
  den1: denominador del primer sumando
  den2: denominador del segundo sumando
  den3: denominador del tercer sumando

Output = aceleracion en una de las componentes para el cuerpo 0
*/
FLOAT ecuacion_gravedad(FLOAT m_1, FLOAT m_2, FLOAT m_3, FLOAT p0, FLOAT p1, FLOAT p2, FLOAT p3, FLOAT den1, FLOAT den2, FLOAT den3);


FLOAT ecuacion_velocidad(FLOAT v_cuerpo);


//Implementacion de Runge_Kutta4 para el un time-step pequeno en una dimension de un solo cuerpo
FLOAT *rungekutta4(FLOAT p_cuerpo, FLOAT p1, FLOAT p2, FLOAT p3,  FLOAT den1, FLOAT den2, FLOAT den3, FLOAT v_cuerpo, FLOAT m1, FLOAT m2, FLOAT m3, FLOAT h);



int main(int argc, char **argv){

  char *filename;
  FLOAT *datos, *dens;
  FLOAT x_sol, y_sol, z_sol, x_tierra, y_tierra, z_tierra, x_luna, y_luna, z_luna, x_ast, y_ast, z_ast;
  FLOAT v_x_sol, v_y_sol, v_z_sol, v_x_tierra, v_y_tierra, v_z_tierra, v_x_luna, v_y_luna, v_z_luna, v_x_ast, v_y_ast, v_z_ast;
  FLOAT masa_sol, masa_luna, masa_tierra, masa_ast;
  FLOAT timestep, year, h, N;
  int i, j, ciclos;
  FLOAT *cx_sol, *cy_sol, *cz_sol, *cx_tierra, *cy_tierra, *cz_tierra, *cx_luna, *cy_luna, *cz_luna, *cx_ast, *cy_ast, *cz_ast;
  //FLOAT d1, d2, d3;

  //Memory allocation de los punteros
  filename = malloc(10*sizeof(char));
  dens = malloc(3*sizeof(FLOAT));
  cx_sol = malloc(2*sizeof(FLOAT));
  cy_sol = malloc(2*sizeof(FLOAT));
  cz_sol = malloc(2*sizeof(FLOAT));
  cx_tierra = malloc(2*sizeof(FLOAT));
  cy_tierra = malloc(2*sizeof(FLOAT));
  cz_tierra = malloc(2*sizeof(FLOAT));
  cx_luna = malloc(2*sizeof(FLOAT));
  cy_luna = malloc(2*sizeof(FLOAT));
  cz_luna = malloc(2*sizeof(FLOAT));
  cx_ast = malloc(2*sizeof(FLOAT));
  cy_ast = malloc(2*sizeof(FLOAT));
  cz_ast = malloc(2*sizeof(FLOAT));

  //Inicializacion de variables que vienen por parametro
  filename = argv[1];
  timestep = atof(argv[2]);
  //year = argv[3];

  h = timestep/4.0; //Como se va a mover cada cuerpo por aparte, el intervalo de movimiento para cada cuerpo sera esto

  //Un anio tiene 31536000 segundos
  N = 31536000/timestep;
  ciclos = (int)100/timestep;

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

  for(i = 0; i<(ciclos); i++){

    //Imprimir archivo con coordenadas  
    FILE *out;
    out = fopen("orbitas.txt", "a");
    if(!out){
      printf("Problemas escribiendo archivo\n");
      exit(1);
    }

    fprintf(out, "%g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g \t %g\n", x_sol, y_sol, z_sol, x_tierra, y_tierra, z_tierra, x_luna, y_luna, z_luna, x_ast, y_ast, z_ast);
   
    fclose(out);

    //Calculo de movimiento del sol
    dens = calcular_distancia(x_sol, y_sol, z_sol, x_tierra, y_tierra, z_tierra, x_luna, y_luna, z_luna, x_ast, y_ast, z_ast);
    cx_sol = rungekutta4(x_sol, x_tierra, x_luna, x_ast, dens[0], dens[1], dens[2], v_x_sol, masa_tierra, masa_luna, masa_ast, h);  
    cy_sol = rungekutta4(y_sol, y_tierra, y_luna, y_ast, dens[0], dens[1], dens[2], v_y_sol, masa_tierra, masa_luna, masa_ast, h);
    cz_sol = rungekutta4(z_sol, z_tierra, z_luna, z_ast, dens[0], dens[1], dens[2], v_z_sol, masa_tierra, masa_luna, masa_ast, h);

    //Calculo del movimiento de la tierra
    dens = calcular_distancia(x_tierra, y_tierra, z_tierra, x_sol, y_sol, z_sol, x_luna, y_luna, z_luna, x_ast, y_ast, z_ast);
    cx_tierra = rungekutta4(x_tierra, x_sol, x_luna, x_ast, dens[0], dens[1], dens[2], v_x_tierra, masa_sol, masa_luna, masa_ast, h);
    cy_tierra = rungekutta4(y_tierra, y_sol, y_luna, y_ast, dens[0], dens[1], dens[2], v_y_tierra, masa_sol, masa_luna, masa_ast, h);
    cz_tierra = rungekutta4(z_tierra, z_sol, z_luna, z_ast, dens[0], dens[1], dens[2], v_z_tierra, masa_sol, masa_luna, masa_ast, h);

    //Calculo del movimiento de la luna
    dens = calcular_distancia(x_luna, y_luna, z_luna, x_sol, y_sol, z_sol, x_tierra, y_tierra, z_tierra, x_ast, y_ast, z_ast);
    cx_luna = rungekutta4(x_luna, x_sol, x_tierra, x_ast, dens[0], dens[1], dens[2], v_x_luna, masa_sol, masa_tierra, masa_ast, h);
    cy_luna = rungekutta4(y_luna, y_sol, y_tierra, y_ast, dens[0], dens[1], dens[2], v_y_luna, masa_sol, masa_tierra, masa_ast, h);
    cz_luna = rungekutta4(z_luna, z_sol, z_tierra, z_ast, dens[0], dens[1], dens[2], v_z_luna, masa_sol, masa_tierra, masa_ast, h);

    //Calculo del movimiento del asteroide
    dens = calcular_distancia(x_ast, y_ast, z_ast, x_sol, y_sol, z_sol, x_tierra, y_tierra, z_tierra, x_luna, y_luna, z_luna);
    cx_ast = rungekutta4(x_ast, x_sol, x_tierra, x_luna, dens[0], dens[1], dens[2], v_x_ast, masa_sol, masa_tierra, masa_luna, h);
    cy_ast = rungekutta4(y_ast, y_sol, y_tierra, y_luna, dens[0], dens[1], dens[2], v_y_ast, masa_sol, masa_tierra, masa_luna, h);
    cz_ast = rungekutta4(z_ast, z_sol, z_tierra, z_luna, dens[0], dens[1], dens[2], v_z_ast, masa_sol, masa_tierra, masa_luna, h);

    //Actualizacion de posiciones y velocidades
    x_sol = cx_sol[0];
    y_sol = cy_sol[0];
    z_sol = cz_sol[0];
    x_tierra = cx_tierra[0];
    y_tierra = cy_tierra[0];
    z_tierra = cz_tierra[0];
    x_luna = cx_luna[0];
    y_luna = cy_luna[0];
    z_luna = cz_luna[0];
    x_ast = cx_ast[0];
    y_ast = cy_ast[0];
    z_ast = cz_ast[0];
    v_x_sol = cx_sol[1];
    v_y_sol = cy_sol[1];
    v_z_sol = cz_sol[1];
    v_x_tierra = cx_tierra[1];
    v_y_tierra = cy_tierra[1];
    v_z_tierra = cz_tierra[1];
    v_x_luna = cx_luna[1];
    v_y_luna = cy_luna[1];
    v_z_luna = cz_luna[1];
    v_x_ast = cx_ast[1];
    v_y_ast = cy_ast[1];
    v_z_ast = cz_ast[1];   

    printf("%d\n", i);

  }

  return 0;

}

FLOAT ecuacion_gravedad(FLOAT m_1, FLOAT m_2, FLOAT m_3, FLOAT p0, FLOAT p1, FLOAT p2, FLOAT p3, FLOAT r1, FLOAT r2, FLOAT r3){

  FLOAT aceleracion;
  //FLOAT G; //Constante de gravitacion universal

  //G = 39.486*31536000;

  //Calculo para componente en deseada0.
  aceleracion = -G_GRAV*((m_1*(p0-p1)/pow(r1, 1.5)) + (m_2*(p0-p2)/pow(r2, 1.5)) + (m_3*(p0-p3)/pow(r3, 1.5)));

  printf("%g\n", aceleracion);

  return aceleracion;
}

FLOAT *calcular_distancia(FLOAT x0, FLOAT y0, FLOAT z0, FLOAT x1, FLOAT y1, FLOAT z1, FLOAT x2, FLOAT y2, FLOAT z2, FLOAT x3, FLOAT y3, FLOAT z3){

  FLOAT *dis;

  dis = malloc(3*sizeof(FLOAT));
  dis[0] = sqrt(pow((x0 - x1), 2 ) + pow ( (y0 - y1), 2 ) + pow((z0 - z1) , 2)); 
  dis[1] = sqrt(pow((x0 - x2), 2 ) + pow ( (y0 - y2), 2 ) + pow((z0 - z2) , 2));
  dis[2] = sqrt(pow((x0 - x3), 2 ) + pow ( (y0 - y3), 2 ) + pow((z0 - z3) , 2));

  return dis;
}

FLOAT *rungekutta4(FLOAT p_cuerpo, FLOAT p1, FLOAT p2, FLOAT p3,  FLOAT den1, FLOAT den2, FLOAT den3, FLOAT v_cuerpo, FLOAT m1, FLOAT m2, FLOAT m3, FLOAT h){

  FLOAT *nuevas_coordenadas;
  FLOAT k1_posicion, k1_velocidad, k2_posicion, k2_velocidad,k3_posicion, k3_velocidad,k4_posicion, k4_velocidad, kprom_posicion, kprom_velocidad;
  FLOAT p_1, v_1, p_2, v_2, p_3, v_3, p_nueva, v_nueva, acel;

  nuevas_coordenadas = malloc(2*sizeof(FLOAT));

  acel = ecuacion_gravedad(m1, m2, m3, p_cuerpo, p1, p2, p3, den1, den2, den3);


  k1_posicion = v_cuerpo;
  k1_velocidad = acel;

  //Primer paso
  p_1 = p_cuerpo + h*0.5*k1_posicion;
  v_1 = v_cuerpo + h*0.5*k1_velocidad;
  k2_posicion = v_1;
  k2_velocidad = acel;

  //Segundo paso
  p_2 = p_cuerpo + h*0.5*k2_posicion;
  v_2 = v_cuerpo + h*0.5*k2_velocidad;
  k3_posicion = v_2;
  k3_velocidad = acel;

  //Tercer paso
  p_3 = p_cuerpo + h*k3_posicion;
  v_3 = v_cuerpo + h*k3_velocidad;
  k4_posicion = v_3;
  k4_velocidad = acel;

  //Cuarto paso
  kprom_posicion = (1.0/6.0)*(k1_posicion + 2.0*k2_posicion + 2.0*k3_posicion + k4_posicion);
  kprom_velocidad = (1.0/6.0)*(k1_velocidad + 2.0*k2_velocidad + 2.0*k3_velocidad + k4_velocidad);

  //Creacion de nuevos valores
  p_nueva = p_cuerpo + h*kprom_posicion;
  v_nueva = v_cuerpo + h*kprom_velocidad;

  nuevas_coordenadas[0] = p_nueva;
  nuevas_coordenadas[1] = v_nueva;

  return nuevas_coordenadas;

} 


FLOAT *importacion_datos(char *filename){

  FILE *in;
  int i;
  char var[50];
  FLOAT *datos;
  char valor[13]; 

  //Memory allocation de los punteros
  datos = malloc(28*sizeof(FLOAT));  

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



  


