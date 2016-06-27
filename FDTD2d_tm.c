/***********************
 *   includes, system   *
 ***********************/
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>

/*******************************
 *   OpenGL Graphics includes   *
 *******************************/
#include <GL/glew.h>
#include <GL/glut.h>


/****************
 *   constants   *
 ****************/
const unsigned int window_x_pos = 128;
const unsigned int window_y_pos = 128;

/* const unsigned int window_width = 256; */
/* const unsigned int window_height = 256; */
const unsigned int window_width = 512;
const unsigned int window_height = 512;


const unsigned int grid_width = 256;
const unsigned int grid_height = 256;
/* const unsigned int grid_width = 512; */
/* const unsigned int grid_height = 512; */


// unsigned int grid_width = 51;
// unsigned int grid_height = 51;

const float min_x = -0.5f;
const float max_x = 0.5f;
const float min_y = -0.5f;
const float max_y = 0.5f;


/****************
 *   parameter   *
 ****************/
float c = 2.99792458e8;
float freq = 1.0e9;
float lambda;
float resolution = 20.0;
float delta_x;
float delta_y;
float alpha = 0.5;
float delta_t;
float step;
float mu0;
float sigma0 = 0;
float epsilon0 = 8.854187e-12f;
float epsilonMax;
int M = 4;
int L = 12;
float r0 = -6;
float ecmax;
float Ez_max = 2.060459378159e-03f;
float Ez_min = -7.196258220476e-04f;
float Ez_range;
float Ez_yellow;
float Ez_green;
float Ez_lightblue;

// float pulse;
float T = 0.0;

float **Ez, **Hx, **Hy, **sigma_M, **epsilon_M, **mu_M;
float *ECX, *ECY;
float **CEZX, **CEZXL, **CHYX, **CHYXL, **CEZY, **CEZYL, **CHXY, **CHXYL;
float **EZX, **EZY, **HXY, **HYX;

float **CEZ, **CEZLX, **CEZLY, **CHXLY, **CHYLX;

// grid_2d
// float *h_F, *h_Fn;
// float grid_dx, grid_dy;

// texture
GLuint tex;
GLubyte *h_g_data;

// counter
unsigned int kt = 1;
unsigned int incrt = 1;

float anim_time = 0.0f;
float anim_dt;

bool flag = false;

#define MAX(a,b) ((a > b) ? a : b)
#define MIN(a,b) ((a < b) ? a : b)

// frame rate
int GLframe = 0;
int GLtimenow = 0;
int GLtimebase = 0;

// mouse controls
int mouse_old_x, mouse_old_y;
int mouse_buttons = 0;
float rotate_x = 0.0f, rotate_y = 0.0f;
// float translate_z = -3.0f;
float translate_z = 0.0f;

// FILE *Ez_op;
//
float r_1=grid_height*1/2, r_2=grid_height*3/4;


/**************************
 *   forward declaration   *
 **************************/

/*** GL functionality ***/
bool initGL(void);

/*** rendering callbacks ***/
void display(void);
// void resize(int width, int height);
void keyboard(unsigned char key, int x, int y);
void mouse(int button, int state, int x, int y);
void motion(int x, int y);
void cleanUp(void);

// void launchCPUKernel(float **Ex, float **Hx, float **sigma_M, float **epsilon_M, float **mu_M, float *ECX, float *ECY, float **CEZX, float **CEZXL, float **CHYX, float **CHYXL, float **CEZY, float **CEZYL, float **CHXY, float **CHXYL, float **EZX, float **EZY, float **HXY, float **HYX, float **CEZ, float **CEZLX, float **CEZLY, float **CHXLY, float **CHYLX, unsigned int width, unsigned int height, float dx, float dy, float kappa, float dt, float max_density, int *g_data);
void launchCPUKernel(void);
void runCPUKernel(void);
// void h_FDTD2d_tm(float *F, float *Fn, unsigned int width, unsigned int height, float c0, float c1, float c2, float max_density, int *g_data);
void h_FDTD2d_tm(void);

void createTexture(GLuint *tex, unsigned int width, unsigned int height);
void deleteTexture(GLuint* tex);

float h_clamp(float x, float a, float b);

void setInitialDensity(unsigned int width, unsigned int height);
void malloc_Initialdata(void);
void free_data(void);

void PEC(void);

void drawCircle(int x0, int y0, int r);
void drawCircle3(int x0, int y0, int r);
void drawCircle4(int x0, int y0, int r);
void drawRact(int r);
void drawRact(int r)
{
  int i, j;
  for(i=0;i<grid_width*3/4+r/2;i++)
  {
    j=grid_height/2-r/2;
    Ez[i][j] = 0.0;
    Hx[i][j] = 0.0;
    Hy[i][j] = 0.0;
  }
  for(i=0;i<grid_width*3/4-r/2;i++)
  {
    j=grid_height/2+r/2;
    Ez[i][j] = 0.0;
    Hx[i][j] = 0.0;
    Hy[i][j] = 0.0;
  }
  for(j=grid_height/2-r/2;j<grid_height;j++)
  {
    i=grid_width*3/4+r/2;
    Ez[i][j] = 0.0;
    Hx[i][j] = 0.0;
    Hy[i][j] = 0.0;
  }
  for(j=grid_height/2+r/2;j<grid_height;j++)
  {
    i=grid_width*3/4-r/2;
    Ez[i][j] = 0.0;
    Hx[i][j] = 0.0;
    Hy[i][j] = 0.0;
  }
}

void drawCircle3(int x0, int y0, int r)
{
  int x=0, y=0;
  for(y=-r; y<=r; y++)
  {
    for(x=-r; x<=r; x++)
    {
      if(x*x+y*y <= r*r)
      {
        if(x0+x>=0 && x0+x<grid_width && y0+y>=0 && y0+y<grid_height)
        {
          Ez[x0+x][y0+y]=0.0;
          Hx[x0+x][y0+y]=0.0;
          Hy[x0+x][y0+y]=0.0;
        }
      }
    }
  }
}
void drawCircle4(int x0, int y0, int r)
{
  int x=0, y=0;
  for(y=-r; y<=r; y++)
  {
    for(x=-r; x<=r; x++)
    {
      if(x*x+y*y > r*r)
      {
        if(x0+x>=0 && x0+x<grid_width && y0+y>=0 && y0+y<grid_height)
        {
          Ez[x0+x][y0+y]=0.0;
          Hx[x0+x][y0+y]=0.0;
          Hy[x0+x][y0+y]=0.0;
        }
      }
    }
  }
}

void drawCircle(int x0, int y0, int r)
{
  int x = r;
  int y = 0;
  int err = 0;
  while(x>=y)
  {
    if(x0+x <= grid_width && y0+y <= grid_height && x0+x >= 0 && y0+y >=0)
      Ez[x0+x][y0+y]=0;
    if(x0+y <= grid_width && y0+x <= grid_height && x0+y >= 0 && y0+x >=0)
      Ez[x0+y][y0+x]=0;
    if(x0-y <= grid_width && y0+x <= grid_height && x0-y >= 0 && y0+x >=0)
      Ez[x0-y][y0+x]=0;
    if(x0-x <= grid_width && y0+y <= grid_height && x0-x >= 0 && y0+y >=0)
      Ez[x0-x][y0+y]=0;
    if(x0-x <= grid_width && y0-y <= grid_height && x0-x >= 0 && y0-y >=0)
      Ez[x0-x][y0-y]=0;
    if(x0-y <= grid_width && y0-x <= grid_height && x0-y >= 0 && y0-x >=0)
      Ez[x0-y][y0-x]=0;
    if(x0+y <= grid_width && y0-x <= grid_height && x0+y >= 0 && y0-x >=0)
      Ez[x0+y][y0-x]=0;
    if(x0+x <= grid_width && y0-y <= grid_height && x0+x >= 0 && y0-y >=0)
      Ez[x0+x][y0-y]=0;

    y+=1;
    err += 1+2*y;
    if(2*(err-x)+1>0)
    {
      x-=1;
      err+=1-2*x;
    }
  }
}

void PEC(void)
{
  /* drawCircle3(0, grid_height, r_1); */
  /* drawCircle4(0, grid_height, r_2); */

  drawRact(30);

}

void h_FDTD2d_tm(void)
{
  unsigned int i, j, index;
  float pulse;
  pulse  =  sin((((kt - 1)%(int)step)+1)*2.0*M_PI/step);

  //Ez
  for(j = 1; j < grid_height-1; j++){
    for(i = 1; i < grid_width-1; i++){
      /* if(i==grid_width/3 && j==grid_height/3){ */
      /* if(i==grid_width/2 && j==grid_height/2){ */
      if(i == grid_width*2/3 && j==grid_height/2){
        Ez[i][j] = 1.0/376.7 * pulse;
      }else{
        Ez[i][j] = CEZ[i][j] * Ez[i][j] + CEZLX[i][j] * (Hy[i][j]-Hy[i-1][j]) - CEZLY[i][j] * (Hx[i][j]-Hx[i][j-1]);
      }
    }
  }
  PEC();

  /* Ez for PML */
  for(j = 1; j<grid_height - 1; j++){
    for(i = 1; i<grid_width - 1; i++){
      if(i<L || (i>grid_width-L-1) || j<L || (j>grid_height-L-1)){
        EZX[i][j]=CEZX[i][j] * EZX[i][j] + CEZXL[i][j] * (Hy[i][j] - Hy[i-1][j]);
        EZY[i][j]=CEZY[i][j] * EZY[i][j] - CEZYL[i][j] * (Hx[i][j] - Hx[i][j-1]);
        Ez[i][j]=EZX[i][j]+EZY[i][j];
      }
    }
  }

  T=T+delta_t/2;

  //Hx
  for(j = 0; j<grid_height - 1; j++){
    for(i = 1; i<grid_width - 1; i++){
      Hx[i][j] = Hx[i][j] - (CHXLY[i][j]*(Ez[i][j+1]-Ez[i][j]));
    }
  }

  /* //Hx for PML*/
  for(j = 0; j<grid_height - 1; j++){
    for(i = 1; i<grid_width - 1; i++){
      if(i<L || i>grid_width-L-1 || j<L || j>grid_height-L-1){
        HXY[i][j]=CHXY[i][j]*HXY[i][j]-CHXYL[i][j]*(Ez[i][j+1]-Ez[i][j]);
        Hx[i][j]=HXY[i][j];
      }
    }
  }

  //Hy
  for(j = 1; j<grid_height - 1; j++){
    for(i = 0; i<grid_width - 1; i++){
      Hy[i][j] = Hy[i][j] + (CHYLX[i][j]*(Ez[i+1][j]-Ez[i][j]));
    }
  }

  //Hy for PML
  for(j = 1; j<grid_height - 1; j++){
    for(i = 0; i<grid_width - 1; i++){
      if(i<L || i>grid_width-L-1 || j<L || j>grid_height-L-1){
        HYX[i][j]=CHYX[i][j]*HYX[i][j]+CHYXL[i][j]*(Ez[i+1][j]-Ez[i][j]);
        Hy[i][j]=HYX[i][j];
      }
    }
  }

  T=T+delta_t/2;

  /***create graphic data***/
  float v;
  for(j=0; j<grid_height; j++){
    for(i=0; i<grid_width; i++){
      index = grid_width * j + i;
      v = Ez[i][j];
      v = h_clamp(v, Ez_min, Ez_max);

      if(v > Ez_yellow) {
        h_g_data[index*3] = (GLubyte)255;
        h_g_data[index*3+1] = (GLubyte)(255-(v-Ez_yellow)/(Ez_max-Ez_yellow)*255);
        h_g_data[index*3+2] = (GLubyte)0;
      }else if(v > Ez_green){
        h_g_data[index*3] = (GLubyte)((v-Ez_green)/(Ez_yellow-Ez_green)*255);
        h_g_data[index*3+1] = (GLubyte)255;
        h_g_data[index*3+2] = (GLubyte)0;
      }else if(v > Ez_lightblue) {
        h_g_data[index*3] = (GLubyte)0;
        h_g_data[index*3+1] = (GLubyte)255;
        h_g_data[index*3+2] = (GLubyte)(255-(v-Ez_lightblue)/(Ez_green-Ez_lightblue)*255);
      }else{
        h_g_data[index*3] = (GLubyte)0;
        h_g_data[index*3+1] = (GLubyte)((v-Ez_min)/(Ez_lightblue-Ez_min)*255);
        h_g_data[index*3+2] = (GLubyte)255;
      }
    }
  }
  kt++;
}

void launchCPUKernel(void)
{
  h_FDTD2d_tm();
}

float h_clamp(float x, float a, float b)
{
  if (x < a)
    x = a;
  if (x > b)
    x = b;
  return x;
}

void malloc_Initialdata(void)
{
  int i;
  Ez  =  (float **)malloc(sizeof(float*) * grid_height);
  Hx  =  (float **)malloc(sizeof(float*) * grid_height);
  Hy  =  (float **)malloc(sizeof(float*) * grid_height);

  sigma_M  =  (float **)malloc(sizeof(float*) * grid_height);
  epsilon_M =  (float **)malloc(sizeof(float*) * grid_height);
  mu_M    =  (float **)malloc(sizeof(float*) * grid_height);

  ECX = (float *)malloc(sizeof(float) * grid_width);
  ECY = (float *)malloc(sizeof(float) * grid_height);

  CEZX  = (float **)malloc(sizeof(float*) * grid_height);
  CEZXL = (float **)malloc(sizeof(float*) * grid_height);
  CHYX  = (float **)malloc(sizeof(float*) * grid_height);
  CHYXL = (float **)malloc(sizeof(float*) * grid_height);
  CEZY  = (float **)malloc(sizeof(float*) * grid_height);
  CEZYL = (float **)malloc(sizeof(float*) * grid_height);
  CHXY  = (float **)malloc(sizeof(float*) * grid_height);
  CHXYL = (float **)malloc(sizeof(float*) * grid_height);

  EZX = (float **)malloc(sizeof(float*) * grid_height);
  EZY  = (float **)malloc(sizeof(float*) * grid_height);
  HXY = (float **)malloc(sizeof(float*) * grid_height);
  HYX = (float **)malloc(sizeof(float*) * grid_height);


  CEZ  = (float **)malloc(sizeof(float*) * grid_height);
  CEZLX = (float **)malloc(sizeof(float*) * grid_height);
  CEZLY  = (float **)malloc(sizeof(float*) * grid_height);
  CHXLY = (float **)malloc(sizeof(float*) * grid_height);
  CHYLX  = (float **)malloc(sizeof(float*) * grid_height);

  for(i=0;i<grid_height;i++){
    Ez[i] = (float *)malloc(sizeof(float) * grid_width);
    Hx[i] = (float *)malloc(sizeof(float) * grid_width);
    Hy[i] = (float *)malloc(sizeof(float) * grid_width);
    sigma_M[i] = (float *)malloc(sizeof(float) * grid_width);
    epsilon_M[i] = (float *)malloc(sizeof(float) * grid_width);
    mu_M[i] = (float *)malloc(sizeof(float) * grid_width);

    CEZX[i]  = (float *)malloc(sizeof(float) * grid_width);
    CEZXL[i] = (float *)malloc(sizeof(float) * grid_width);
    CHYX[i]  = (float *)malloc(sizeof(float) * grid_width);
    CHYXL[i] = (float *)malloc(sizeof(float) * grid_width);
    CEZY[i]  = (float *)malloc(sizeof(float) * grid_width);
    CEZYL[i] = (float *)malloc(sizeof(float) * grid_width);
    CHXY[i]  = (float *)malloc(sizeof(float) * grid_width);
    CHXYL[i] = (float *)malloc(sizeof(float) * grid_width);

    EZX[i] = (float *)malloc(sizeof(float) * grid_width);
    EZY[i] = (float *)malloc(sizeof(float) * grid_width);
    HXY[i] = (float *)malloc(sizeof(float) * grid_width);
    HYX[i] = (float *)malloc(sizeof(float) * grid_width);

    CEZ[i] =  (float *)malloc(sizeof(float) * grid_width);
    CEZLX[i] = (float *)malloc(sizeof(float) * grid_width);
    CEZLY[i] = (float *)malloc(sizeof(float) * grid_width);
    CHXLY[i] = (float *)malloc(sizeof(float) * grid_width);
    CHYLX[i] = (float *)malloc(sizeof(float) * grid_width);
  }
}

void free_data(void)
{
  free(Ez);
  free(Hx);
  free(Hy);
  free(ECX);
  free(ECY);
  free(CEZX);
  free(CEZXL);
  free(CHYX);
  free(CHYXL);
  free(CEZY);
  free(CEZYL);
  free(CHXY);
  free(CHXYL);
  free(EZX);
  free(EZY);
  free(HXY);
  free(HYX);

  free(CEZ);
  free(CEZLX);
  free(CEZLY);
  free(CHXLY);
  free(CHYLX);
}

/************************
 * Initialize Data      *
 ************************/
void setInitialData(unsigned int width, unsigned int height)
  // unsigned int width, height; 格子のX方向とY方向の解像度．
  // float *max_density; 濃度の最大値．
{
  /*unsigned int i, j, index;
    float x, y;
    grid_dx = (max_x - min_x) / width;
    grid_dy = (max_y - min_y) / height;

  // 矩形の濃度場を与える．
   *max_density = 1.0f;
   h_F = (float *)malloc(width * height * sizeof(float));
   h_Fn = (float *)malloc(width * height * sizeof(float));
   for (i = 0; i < height; i++) {
   for (j = 0; j < width; j++) {
   index = width * i + j;
   x = grid_dx * (j + 0.5f) + min_x;
   y = grid_dy * (i + 0.5f) + min_y;
   if ((x > (0.75f * min_x + 0.25f * max_x)) && (x < (0.25f * min_x + 0.75f * max_x)) && (y > (0.75f * min_y + 0.25f * max_y)) && (y < (0.25f * min_y + 0.75f * max_y))){
   h_F[index] = *max_density;
   }else{
   h_F[index] = 0.0f;
   }
   }
   }
   */

  lambda = c / freq;
  delta_x = lambda / resolution;
  delta_y = lambda / resolution;	
  delta_t = (1.0 / (sqrt(pow((1 / delta_x), 2.0)+pow((1 / delta_y), 2.0))))*(1.0 / c)*alpha;
  anim_dt = (1.0 / (sqrt(pow((1 / delta_x), 2.0)+pow((1 / delta_y), 2.0))))*(1.0 / c)*alpha;
  step = 1.0 / freq / delta_t;
  mu0 = 1.0e-7f * 4.0 * M_PI;
  ecmax = -(M+1)*epsilon0*c / (2.0*L*delta_x)*r0;
  Ez_range = Ez_max-Ez_min; // 2.7800852e-03f 
  Ez_yellow = Ez_range*0.75f+Ez_min;
  Ez_green = Ez_range*0.50f+Ez_min;
  Ez_lightblue = Ez_range*0.25f+Ez_min;

  // printf("%.12e\n", Ez_lightblue);

  int i, j;
  float Z, ZZ;
  for(j = 0; j<grid_height; j++){
    for(i = 0; i<grid_width; i++){
      mu_M[i][j]  =  mu0;
      epsilon_M[i][j] = epsilon0;
      sigma_M[i][j] = sigma0;
    }
  }

  for(j = 0; j<grid_height; j++){
    for(i = 0;i<grid_width; i++){
      Ez[i][j] = 0.0;
      Hx[i][j] = 0.0;
      Hy[i][j] = 0.0;
      CEZX[i][j] = 0.0;
      CEZXL[i][j] = 0.0;
      CHYX[i][j] = 0.0;
      CHYXL[i][j] = 0.0;
      CEZY[i][j] = 0.0;
      CEZYL[i][j] = 0.0;
      CHXY[i][j] = 0.0;
      CHXYL[i][j] = 0.0;

      CEZ[i][j] = 0.0;
      CEZLX[i][j]=0.0;
      CEZLY[i][j]=0.0;
      CHXLY[i][j]=0.0;
      CHYLX[i][j]=0.0;
    }
  }

  for(i=0;i<grid_width;i++){
    ECX[i]=0.0;
  }

  for(j=0;j<grid_width;j++){
    ECY[j]=0.0;
  }

  for(i=0;i<L;i++){
    ECX[i] = ecmax * pow((L-i+0.5)/L,M);
    ECX[grid_width-i-1] = ECX[i];
    ECY[i] = ECX[i];
    ECY[grid_height-i-1] = ECX[i];
  }

  //PML init
  for(i=0;i<grid_width;i++){
    for(j=0;j<grid_height;j++){
      Z = (ECX[i] * delta_t)/(2.0*epsilon_M[i][j]);
      CEZX[i][j]=(1-Z)/(1+Z);
      CEZXL[i][j]=(delta_t/epsilon_M[i][j])/(1+Z)*(1.0/delta_x);
      CHYX[i][j]=(1-Z)/(1+Z);
      CHYXL[i][j]=(delta_t/mu_M[i][j])*(1.0/delta_x);
      Z = (ECY[j]*delta_t)/(2.0*epsilon_M[i][j]);
      CEZY[i][j]=(1-Z)/(1+Z);
      CEZYL[i][j]=(delta_t/epsilon_M[i][j])/(1+Z)*(1.0/delta_y);
      CHXY[i][j]=(1-Z)/(1+Z);
      CHXYL[i][j]=(delta_t/mu_M[i][j])*(1.0/delta_y);
    }
  }

  //FDTD init
  for(i=0;i<grid_width;i++){
    for(j=0;j<grid_height;j++){
      ZZ = (sigma_M[i][j] * delta_t)/(2.0*epsilon_M[i][j]);
      CEZ[i][j]=(1-ZZ)/(1+ZZ);
      CEZLX[i][j]=(delta_t/epsilon_M[i][j])/(1+ZZ)*(1.0/delta_x);
      CEZLY[i][j]=(delta_t/epsilon_M[i][j])/(1+ZZ)*(1.0/delta_y);
      CHXLY[i][j]=delta_t/mu_M[i][j]*(1.0/delta_y);
      CHYLX[i][j]=delta_t/mu_M[i][j]*(1.0/delta_x);
    }
  }
}

// Create texture
void createTexture(GLuint *tex, unsigned int width, unsigned int height)
  // GLuint *tex; texture
  // unsinged int width, height; テクスチャに格納するイメージデータの解像度．
{
  glGenTextures(1, tex);

  // 基本的なパラメータの設定．
  glBindTexture(GL_TEXTURE_2D, *tex);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

  // 格納するイメージデータの定義．
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);
}

// Delete texture
void deleteTexture(GLuint* tex)
  // GLuint *tex; テクスチャ．
{
  glDeleteTextures(1, tex);
  *tex = 0;
}

/********************
 *   Initialize GL   *
 ********************/
bool initGL(void){
  // glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
  // glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  glClearColor(135.0f/255.0f, 206.0f/255.0f, 250.0f/255.0f, 1.0f);
  glEnable(GL_TEXTURE_2D);


  /*** initialize glew ***/
  glewInit();
  /*if (!glewIsSupported("GL_VERSION_2_0")) {
    fprintf(stderr, "ERROR: Support for necessary OpenGL extensions missing.");
    return false;
    }*/

  // create texture
  // h_g_data = (int *)malloc(grid_width * grid_height * sizeof(int));
  // h_g_data = new GLubyte[grid_width * grid_height * 3];
  h_g_data = (GLubyte *)malloc(sizeof(GLubyte) * grid_width * grid_height * 3);
  createTexture(&tex, grid_width, grid_height);

  return true;
}

/********************
 *   Run CPUKernel   *
 ********************/
void runCPUKernel(void) 
{
  // CPU処理の呼び出し．

  launchCPUKernel();
  // h_FDTD2d_tm();

  anim_time += anim_dt;

  // 画像データをテクスチャへ転送．	
  glBindTexture(GL_TEXTURE_2D, tex);
  glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, grid_width, grid_height, GL_RGB, GL_UNSIGNED_BYTE, h_g_data);
}

/************************
 *   Display callback    *
 ************************/
void display(void)
{
  /*** update density ***/
  if(flag == true){
    runCPUKernel();
  }

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(- 1.0, 1.0, - 1.0, 1.0, -1.0, 1.0);

  /*** translate and rotate ***/
  glTranslatef(0.0, 0.0, translate_z);
  glRotatef(rotate_x, 1.0, 0.0, 0.0);
  glRotatef(rotate_y, 0.0, 1.0, 0.0);

  /*** draw square ***/
  glClear(GL_COLOR_BUFFER_BIT);
  glBindTexture(GL_TEXTURE_2D, tex);
  glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
  glBegin(GL_QUADS);
  glTexCoord2f(0.0f, 0.0f); glVertex2f(- 1.0, - 1.0); 
  glTexCoord2f(1.0f, 0.0f); glVertex2f(1.0, - 1.0); 
  glTexCoord2f(1.0f, 1.0f); glVertex2f(1.0, 1.0); 
  glTexCoord2f(0.0f, 1.0f); glVertex2f(- 1.0, 1.0);
  // glTexCoord2f(0.0f, 0.0f); glVertex2f(- 0.5, - 0.5); 
  // glTexCoord2f(0.5f, 0.0f); glVertex2f(0.5, - 0.5); 
  // glTexCoord2f(0.5f, 0.5f); glVertex2f(0.5, 0.5); 
  //glTexCoord2f(0.0f, 0.5f); glVertex2f(- 0.5, 0.5);
  glEnd();	

  /*** update graphics ***/
  glutSwapBuffers();
  glutPostRedisplay();

  if(flag == true){
    incrt++;
    if (incrt % 100 == 0){
      printf("time(%4d): %7.5f\n", incrt, anim_time);
    }
  }
  if(incrt == 400){
    flag = false;
  }

  /*** timer ***/
  GLframe++;
  GLtimenow = glutGet(GLUT_ELAPSED_TIME);
  if(GLtimenow - GLtimebase > 1000)
  {
    char fps[256];
    sprintf(fps, "FDTD2d TM-mode (PBO): %.2f fps", GLframe*1000.0/(GLtimenow-GLtimebase));
    glutSetWindowTitle(fps);

    GLtimebase = GLtimenow;
    GLframe = 0;
  }
}

/****************************
 *   Keyboard events handler *
 ****************************/
void keyboard(unsigned char key, int x, int y)
{
  switch (key){
    case 'x':
      r_1-=0.5;
      r_2+=0.5;
      break;
    case 'z':
      r_1+=0.5;
      r_2-=0.5;
      break;
    case (27):
      exit(EXIT_SUCCESS);
      break;
    case 's':
      flag = !flag;
      break;
    case 'i':
      rotate_x +=2;
      break;
    case 'k':
      rotate_x -=2;
      break;
    case 'j':
      rotate_y +=2;
      break;
    case 'l':
      rotate_y -=2;
      break;
    case 'u':
      translate_z +=0.1;
      break;
    case 'o':
      translate_z -=0.1;
      break;
    case 'I':
      rotate_x +=0.4;
      break;
    case 'K':
      rotate_x -=0.4;
      break;
    case 'J':
      rotate_y +=0.4;
      break;
    case 'L':
      rotate_y -=0.4;
      break;
    case 'U':
      translate_z +=0.02;
      break;
    case 'O':
      translate_z -=0.02;
      break;
    case 'r':
      flag = false;
      incrt = 1;
      kt = 1;
      T = 0.0;
      anim_time = 0.0f;
      rotate_x=0.0f;
      rotate_y=0.0f;
      // translate_z = -3.0f;
      translate_z = 0.0f;
      setInitialData(grid_width, grid_height);
      // printf("Max: %f\n", max_density);
      glutSwapBuffers();
      glutPostRedisplay();
      break;
  }
}

/****************************
 *   Mouse events handler *
 ****************************/
void mouse(int button, int state, int x, int y)
{
  if (state == GLUT_DOWN)
  {
    mouse_buttons |= 1<<button;
  }
  else if (state == GLUT_UP)
  {
    mouse_buttons = 0;
  }

  mouse_old_x = x;
  mouse_old_y = y;
}

void motion(int x, int y)
{
  float dx, dy;
  dx = (float)(x - mouse_old_x);
  dy = (float)(y - mouse_old_y);

  if (mouse_buttons & 1)
  {
    rotate_x += dy * 0.2f;
    rotate_y += dx * 0.2f;
  }
  else if (mouse_buttons & 4)
  {
    translate_z += dy* 0.01f;
  }

  mouse_old_x = x;
  mouse_old_y = y;
}

void cleanUp(void)
{
  // free(h_F);
  // free(h_Fn);

  /*** delete texture ***/
  // fclose(Ez_op);
  deleteTexture(&tex);
  free(h_g_data);
  free_data();
}

/*******************
 *   Program Main   *
 *******************/
int main(int argc, char **argv)
{
  /*
     if((Ez_op=fopen("Ez.txt","a"))==NULL){
     printf("file error\n");
     return 1;
     }
     */
  // printf("%d\n", grid_width);

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
  glutInitWindowPosition(window_x_pos, window_y_pos);
  glutInitWindowSize(window_width, window_height);
  glutCreateWindow("FDTD2d TM-mode (PBO)");
  glutDisplayFunc(display);
  // glutReshapeFunc(resize);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  glViewport(0, 0, window_width, window_height);
  atexit(cleanUp);

  /*** init data ***/
  malloc_Initialdata();
  setInitialData(grid_width, grid_height);
  // printf("Max: %f\n", max_density);

  // anim_dt = 0.2 * MIN(grid_dx*grid_dx, grid_dy*grid_dy) / kappa;

  if(!initGL()){
    return 1;
  }

  glutMainLoop();
  return 0;
}
