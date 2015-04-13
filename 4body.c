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
  
  if(!in)
    {
      printf("Problemas leyendo archivo %s \n", filename);
      exit(1);
      
    }

  for(i=0;i<21;i++)
    {
      fscanf(in, "%s\n", %var);
      printf()

    }

  char buf[50];

  file = fopen(filename,"r");

  while 
  

//}

