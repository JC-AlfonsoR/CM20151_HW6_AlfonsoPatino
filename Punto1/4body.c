#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//Funcion que importa los datos del archivo con nombre que se le pasa como parametro.
double * importacion_datos(char *filename);
/*
La funcion calcula los vectores de posicion inciales y la velocidad inicial para los cuatro cuerpos. Se debe poner como input los datod
importados del archivo que contiene las condiciones iniciales.
*/
double * IC(double *datos);
//double * rungekutta4(double y_0, double *t, double (*f)(*double),  ); 

int main(int argc, char **argv){

  char *filename;
  double *datos;
  double *condiciones_iniciales;

  //Memory allocation de los punteros
  filename = malloc(10*sizeof(char));

  //Inicializacion de variables
  filename = argv[1];


  //Se llama la funcion para importar los datos
  datos  = importacion_datos(filename);

  //Se crean las condiciones iniciales usando los datos importados
  condiciones_iniciales = IC(datos);


  return 0;

}


double * importacion_datos(char *filename){

  FILE *in;
  int i;
  char var[50];
  double *datos;
  char valor[10]; 

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

double * IC(double *datos)
{
  double *IC;
  double r_0_sol;
  double r_0_tierra;
  double r_0_luna;
  double r_0_ast;
  double v_0_sol;
  double v_0_tierra;
  double v_0_luna;
  double v_0_ast;

  //Memory allocation del vector de condiciones iniciales
  IC = malloc(8*sizeof(double));

  //Condiciones iniciales de posiciones
  r_0_sol = pow((pow(datos[0], 2) + pow(datos[1], 2) + pow(datos[2], 2)), 0.5);
  IC[0] = r_0_sol;

  r_0_tierra = pow((pow(datos[3], 2) + pow(datos[4], 2) + pow(datos[5], 2)), 0.5);
  IC[1] = r_0_tierra;

  r_0_luna = pow((pow(datos[6], 2) + pow(datos[7], 2) + pow(datos[8], 2)), 0.5);
  IC[2] = r_0_luna;

  r_0_ast = pow((pow(datos[9], 2) + pow(datos[10], 2) + pow(datos[11], 2)), 0.5);
  IC[3] = r_0_ast;

  //Condiciones iniciales de velocidad
  v_0_sol = pow((pow(datos[12], 2) + pow(datos[13], 2) + pow(datos[14], 2)), 0.5);
  IC[4] = r_0_sol;

  v_0_tierra = pow((pow(datos[15], 2) + pow(datos[16], 2) + pow(datos[17], 2)), 0.5);
  IC[5] = r_0_sol;

  v_0_luna = pow((pow(datos[18], 2) + pow(datos[19], 2) + pow(datos[20], 2)), 0.5);
  IC[6] = r_0_luna;

  v_0_ast = pow((pow(datos[21], 2) + pow(datos[22], 2) + pow(datos[23], 2)), 0.5);
  IC[7] = r_0_ast;


  return IC;
}

  


