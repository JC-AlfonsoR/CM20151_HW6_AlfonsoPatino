#include <stdlib.h>
#include <stdio.h>

void importar_datos(char *filename);

int main(int argc, char **argv){

  char *filename;

  //Memory allocation de los punteros
  filename = malloc(10*sizeof(char));

  //Inicializacion de variables
  filename = argv[1];


  //Se llama la funcion para importar los datos
  importar_datos(filename);

  return 0;

}


importar_datos(char *filename){

  FILE *in;
  int i;
  char var[50];
  double *datos;

  datos = malloc(20*sizeof(double));

  in = fopen(filename, "r");
  
  if(!in)
    {
      printf("Problemas leyendo archivo %s \n", filename);
      exit(1);
      
    }

  for(i=0;i<42;i++)
    {
      fscanf(in, "%s \n", &var);
      printf("%s \n", var);

    }
  fclose(in);
}
  


