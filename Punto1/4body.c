#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define G_GRAV 6.67384E-11;

/*
Se adopta la siguiente convencion de indices. 0: es para el cuerpo que se mueve. 1,2,3 son para los cuerpos que estan afectando al cuerpo
en movimiento
*/

/*
Funcion que importa los datos desde el archivo de condiciones iniciales.
Input: nombre del archivo
Output: array con condiciones iniciales. (Para ver a que corresponde cada valor en el array ver archivo ic.txt)
*/
double * importacion_datos(char *filename);

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
 double *calcular_distancia(double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3);

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
double ecuacion_gravedad(double m_1, double m_2, double m_3, double p0, double p1, double p2, double p3, double den1, double den2, double den3);


double ecuacion_velocidad(double v_cuerpo);


//Implementacion de Runge_Kutta4 para el un time-step pequeno en una dimension de un solo cuerpo
double *rungekutta4(double p_cuerpo, double p1, double p2, double p3,  double den1, double den2, double den3, double v_cuerpo, double m1, double m2, double m3, double h);



int main(int argc, char **argv){

  char *filename;
  double *datos, *dens;
  double x_sol, y_sol, z_sol, x_tierra, y_tierra, z_tierra, x_luna, y_luna, z_luna, x_ast, y_ast, z_ast;
  double v_x_sol, v_y_sol, v_z_sol, v_x_tierra, v_y_tierra, v_z_tierra, v_x_luna, v_y_luna, v_z_luna, v_x_ast, v_y_ast, v_z_ast;
  double masa_sol, masa_luna, masa_tierra, masa_ast;
  double timestep, year, h, N;
  int i, j, ciclos;
  double *cx_sol, *cy_sol, *cz_sol, *cx_tierra, *cy_tierra, *cz_tierra, *cx_luna, *cy_luna, *cz_luna, *cx_ast, *cy_ast, *cz_ast;
  //double d1, d2, d3;

  //Memory allocation de los punteros
  filename = malloc(10*sizeof(char));
  dens = malloc(3*sizeof(double));
  cx_sol = malloc(2*sizeof(double));
  cy_sol = malloc(2*sizeof(double));
  cz_sol = malloc(2*sizeof(double));
  cx_tierra = malloc(2*sizeof(double));
  cy_tierra = malloc(2*sizeof(double));
  cz_tierra = malloc(2*sizeof(double));
  cx_luna = malloc(2*sizeof(double));
  cy_luna = malloc(2*sizeof(double));
  cz_luna = malloc(2*sizeof(double));
  cx_ast = malloc(2*sizeof(double));
  cy_ast = malloc(2*sizeof(double));
  cz_ast = malloc(2*sizeof(double));

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
    printf("%g, %g, %g, %g\n", x_sol, x_tierra, x_luna, x_ast);
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

double ecuacion_gravedad(double m_1, double m_2, double m_3, double p0, double p1, double p2, double p3, double r1, double r2, double r3){

  double aceleracion;
  double G; //Constante de gravitacion universal

  //G = 6.67384E-11;
  G = 39.486;

  //Calculo para componente en deseada
  aceleracion = -G*((m_1*(p0-p1)/pow(r1, 1.5)) + (m_2*(p0-p2)/pow(r2, 1.5)) + (m_3*(p0-p3)/pow(r3, 1.5)));

  printf("%g\n", aceleracion);

  return aceleracion;
}

double *calcular_distancia(double x0, double y0, double z0, double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3){

  double *dis;

  dis = malloc(3*sizeof(double));
  dis[0] = sqrt(pow((x0 - x1), 2 ) + pow ( (y0 - y1), 2 ) + pow((z0 - z1) , 2)); 
  dis[1] = sqrt(pow((x0 - x2), 2 ) + pow ( (y0 - y2), 2 ) + pow((z0 - z2) , 2));
  dis[2] = sqrt(pow((x0 - x3), 2 ) + pow ( (y0 - y3), 2 ) + pow((z0 - z3) , 2));

  return dis;
}

double *rungekutta4(double p_cuerpo, double p1, double p2, double p3,  double den1, double den2, double den3, double v_cuerpo, double m1, double m2, double m3, double h){

  double *nuevas_coordenadas;
  double k1_posicion, k1_velocidad, k2_posicion, k2_velocidad,k3_posicion, k3_velocidad,k4_posicion, k4_velocidad, kprom_posicion, kprom_velocidad;
  double p_1, v_1, p_2, v_2, p_3, v_3, p_nueva, v_nueva;

  nuevas_coordenadas = malloc(2*sizeof(double));

  k1_posicion = v_cuerpo;
  k1_velocidad = ecuacion_gravedad(m1, m2, m3, p_cuerpo, p1, p2, p3, den1, den2, den3);

  //Primer paso
  p_1 = p_cuerpo + h*0.5*k1_posicion;
  v_1 = v_cuerpo + h*0.5*k1_velocidad;
  k2_posicion = v_1;
  k2_velocidad = ecuacion_gravedad(m1, m2, m3, p_1, p1, p2, p3, den1, den2, den3);

  //Segundo paso
  p_2 = p_cuerpo + h*0.5*k2_posicion;
  v_2 = v_cuerpo + h*0.5*k2_velocidad;
  k3_posicion = v_2;
  k3_velocidad = ecuacion_gravedad(m1, m2, m3, p_2, p1, p2, p3, den1, den2, den3);

  //Tercer paso
  p_3 = p_cuerpo + h*k3_posicion;
  v_3 = v_cuerpo + h*k3_velocidad;
  k4_posicion = v_3;
  k4_velocidad = ecuacion_gravedad(m1, m2, m3, p_3, p1, p2, p3, den1, den2, den3);

  //Cuarto paso
  kprom_posicion = (1.0/6.0)*(k1_posicion + 2.0*k2_posicion + 2.0*k3_posicion + k4_posicion);
  kprom_velocidad = (1.0/6.0)*(k1_velocidad + 2.0*k2_velocidad + 2.0*k3_velocidad + k4_velocidad);

  //Creacion de nuevos valores
  p_nueva = p_cuerpo + h*kprom_posicion;
  v_nueva = v_cuerpo + v_cuerpo + h*kprom_velocidad;

  nuevas_coordenadas[0] = p_nueva;
  nuevas_coordenadas[1] = v_nueva;

  return nuevas_coordenadas;

} 


double *importacion_datos(char *filename){

  FILE *in;
  int i;
  char var[50];
  double *datos;
  char valor[13]; 

  //Memory allocation de los punteros
  datos = malloc(28*sizeof(double));  

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



  


