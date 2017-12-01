// prototypes
#include "real.h"

__global__
void bcTB(real *u,real *v,real *w, real *t2, real *t1);

__global__
void bcEW(real *u,real *v,real *w,real *t2, real *t1);

__global__
void bcNS(real *u,real *v,real *w,real *t2, real *t1);

__global__
void pgf2(real *p3,const real *u,const real *v,const real *w,const real *p1,const real *ro_u,const real *ro_w) ;


__global__
void pgf1(real *u,real *v,real *w,const real *t,const real *p1,const real *ro_u,const real *ro_w);

__global__
void diffusion(real *t,real *u,real *v,real *w,const real *tt,const real *uu,const real *vv,const real *ww);

__global__
void advecUVW(real *u,real *v,real *w, const real *uu,const real *vv,const  real *ww);

__global__
void advectionTB(real *tTemp,const real *t2,const real *w);

void update(real **cube1, real **cube2);


__global__
void advectionNS(real *tTemp,const real *t2,const real *v);

__global__
void advectionEW(real *t2,const real *t1,const real *u);


//void advectionEW(real ***q2,real ***u, real ***q1, real dy,real dt,int i1,int i2,int j1,int j2,int k1,int k2,int nydim,real tstep, char advection_type);


void ic(real ***t ,real ***u, real ***v,real ***w ,real *ro_u,real *ro_w,real dx,real dy,real  dz, real deltaU, int i1,int i2,int j1,int j2,int k1,int k2,int bcW,int nx,int ny,int nz,real *x0,real *y0,real *z0, real *deltaTheta,real *deltaV,real thetaBar, real *rx, real *ry, real *rz,real g);



void putfield(char *name,float datatime,float *field,int nx,int ny,int nz);

//void stats(float **q,int i1,int i2,int j1,int j2,int n,float *qmax, float *qmin);
//void plot1d(float *q,int i1,int i2,int nx,int n,float qmax,int overlay,float *qtrue,char *title,char *name);
//void sfc(float **q,int nx,int ny,float time,float angh,float angv,char *label,char *name);
//void contr(float **z,int nx,int ny,float cint,float time,char *label,int colors,int pltzero,int nestX1,int nestX2,int nestY1,int nestY2,char *name);
//void getTakacsError(float **q1,float dx,float dy,int i1,int i2,int j1,int j2,int nx,int ny,float x0,float y0,float r);




// end  of prototypes

