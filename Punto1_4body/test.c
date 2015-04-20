#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define FLOAT float
#define PI 3.141592653589793
#define G_GRAV 39.486//31536000 //units of ua+3 msun-1 s-1

FLOAT * get_memory(int n_points);
void initialize_pos(FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *ic);
void initialize_vel(FLOAT *vx, FLOAT *vy, FLOAT *vz, FLOAT *ic);
void initialize_mass(FLOAT *mass, FLOAT *ic);
void print_status(FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *vx, FLOAT *vy, FLOAT *vz, FLOAT *ax, FLOAT *ay, FLOAT *az, int n_points);
void get_acceleration(FLOAT *ax, FLOAT *ay, FLOAT *az, FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *mass, int n_points);
FLOAT *rungekutta4(FLOAT p_cuerpo, FLOAT vel, FLOAT acel, FLOAT h);
FLOAT * importacion_datos(char *filename);

int main(int argc, char **argv){
  
  /*positions of all particles*/
  FLOAT *x;
  FLOAT *y;
  FLOAT *z;
  
  /*velocities of all particles*/
  FLOAT *v_x;
  FLOAT *v_y;
  FLOAT *v_z;

  /*accelerations of all particles*/
  FLOAT *a_x;
  FLOAT *a_y;
  FLOAT *a_z;

  /*masses*/
  FLOAT *mass;

  /*timestep variables*/
  FLOAT delta_t= 0.001;
  int n_steps = (int)(100.0/delta_t);
  int n_points = 4;
  FLOAT radius = 100.0;
  FLOAT unit_mass = 1.0; 
  FLOAT *ic;
  int i,j;
  
  /*memory allocation*/
  x = get_memory(n_points);
  y = get_memory(n_points);
  z = get_memory(n_points);
  v_x = get_memory(n_points);
  v_y = get_memory(n_points);
  v_z = get_memory(n_points);
  a_x = get_memory(n_points);
  a_y = get_memory(n_points);
  a_z = get_memory(n_points);
  mass = get_memory(n_points);
  ic = get_memory(28);

  ic =  importacion_datos(argv[1]);

  initialize_pos(x,y,z, ic);
  initialize_vel(v_x,v_y,v_z, ic);
  initialize_mass(mass, ic);

  /*implementation of a simple Euler integration*/
  for(i=0;i<n_steps;i++){
    get_acceleration(a_x, a_y, a_z, x, y, z, mass, n_points);
    for(j=0;j<n_points;j++){
      FLOAT *ks;
      ks = malloc(2*sizeof(FLOAT));
      ks = rungekutta4(x[j] , v_x[j], a_x[j], delta_t);
      x[j] = x[j] + delta_t * ks[0];
      v_x[j] = v_x[j] + delta_t * ks[1];

      ks = rungekutta4(y[j], v_y[j], a_y[j], delta_t);
      y[j] = y[j] + delta_t * ks[0];
      v_y[j] = v_y[j] + delta_t * ks[1];

      ks = rungekutta4(z[j], v_z[j], a_z[j], delta_t);
      z[j] = z[j] + delta_t * ks[0];   
      v_z[j] = v_z[j] + delta_t * ks[1] ;

    }
    print_status(x,y,z,v_x,v_y,v_z, a_x, a_y, a_z, n_points);
  }  
}

void get_acceleration(FLOAT *ax, FLOAT *ay, FLOAT *az, FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *mass, int n_points){
  int i,j;
  FLOAT r_ij;
  for(i=0;i<n_points;i++){
    ax[i]=0.0;
    ay[i]=0.0;
    az[i]=0.0;
    
    for(j=0;j<n_points;j++){
      if(j!=i){
	r_ij = (pow((x[i] - x[j]),2.0) +
		pow((y[i] - y[j]),2.0) +
		pow((z[i] - z[j]),2.0));
	r_ij = sqrt(r_ij);
 
	ax[i] += -G_GRAV *mass[j]/ pow(r_ij,1.5) * (x[i] - x[j]);
	ay[i] += -G_GRAV *mass[j]/ pow(r_ij,1.5) * (y[i] - y[j]);
	az[i] += -G_GRAV *mass[j] / pow(r_ij,1.5) * (z[i] - z[j]);
      }
    }
  }  
}

void initialize_pos(FLOAT *x, FLOAT *y, FLOAT *z,  FLOAT  *ic
  ){

  x[0] = ic[0];
  y[0] = ic[1];
  z[0] = ic[2];
  x[1] = ic[3];
  y[1] = ic[4];
  z[1] = ic[5];
  x[2] = ic[6];
  y[2] = ic[7];
  z[2] = ic[8];
  x[3] = ic[9];
  y[3] = ic[10];
  z[3] = ic[1];
}

void initialize_vel(FLOAT *vx, FLOAT *vy, FLOAT *vz, FLOAT *ic){

  vx[0] = ic[12];
  vy[0] = ic[13];
  vz[0] = ic[14];
  vx[1] = ic[15];
  vy[1] = ic[16];
  vz[1] = ic[17];
  vx[2] = ic[18];
  vy[2] = ic[19];
  vz[2] = ic[20];
  vx[3] = ic[21];
  vy[3] = ic[22];
  vz[3] = ic[23];

    

}

void initialize_mass(FLOAT *mass, FLOAT *ic){
  
  mass[0] = ic[24];
  mass[1] = ic[25];
  mass[2] = ic[26];
  mass[3] = ic[27];
}

FLOAT * get_memory(int n_points){
  FLOAT * x; 
  if(!(x = malloc(sizeof(FLOAT) * n_points))){
    printf("problem with memory allocation");
    exit(1);
  }
  return x;
}

void print_status(FLOAT *x, FLOAT *y, FLOAT *z, FLOAT *vx, FLOAT *vy, FLOAT *vz, FLOAT *ax, FLOAT *ay, FLOAT *az, int n_points){
  
  int i;
  FILE *out;
  out = fopen("orbitas.txt", "a");
  for(i=0;i<n_points;i++){
    fprintf(out, "%f %f %f ",
     x[i], y[i], z[i]);
  }
  fprintf(out, "\n");
  fclose(out);
}

FLOAT *rungekutta4(FLOAT p_cuerpo, FLOAT vel, FLOAT acel, FLOAT h){

  FLOAT k1_posicion, k1_velocidad, k2_posicion, k2_velocidad,k3_posicion, k3_velocidad,k4_posicion, k4_velocidad, kprom_posicion, kprom_velocidad;
  FLOAT p_1, v_1, p_2, v_2, p_3, v_3;
  FLOAT *ks;

  ks = malloc(2*sizeof(FLOAT));


  k1_posicion = vel;
  k1_velocidad = acel;

  //Primer paso
  p_1 = p_cuerpo + h*0.5*k1_posicion;
  v_1 = vel + h*0.5*k1_velocidad;
  k2_posicion = v_1;
  k2_velocidad = acel;
  //Segundo paso
  p_2 = p_cuerpo + h*0.5*k2_posicion;
  v_2 = vel + h*0.5*k2_velocidad;
  k3_posicion = v_2;
  k3_velocidad = acel;
  //Tercer paso
  p_3 = p_cuerpo + h*k3_posicion;
  v_3 = vel + h*k3_velocidad;
  k4_posicion = v_3;
  k4_velocidad = acel;

  //Cuarto paso
  ks[0] = (1.0/6.0)*(k1_posicion + 2.0*k2_posicion + 2.0*k3_posicion + k4_posicion);
  ks[1] = (1.0/6.0)*(k1_velocidad + 2.0*k2_velocidad + 2.0*k3_velocidad + k4_velocidad);

  return ks;
}

FLOAT *importacion_datos(char *filename){

  FILE *in;
  int i;
  char var[50];
  FLOAT *datos;
  char valor[13]; 

  //Memory allocation de los punteros
  datos = get_memory(28);  

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