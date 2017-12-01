// prototypes
#include "real.h"

void pgf1(real *** u,real *** v,real *** w,real *** t,real *** p1,real *ro_u,real *ro_w,real dx,real dy,real dz,int i1,int i2,int j1,int j2,int k1,int k2,real tstep,real g,real thetaBar, int bcw) ;

void pgf2(real *** p3,real *** u,real *** v,real *** w,real *** t,real *** p1,real *ro_u,real *ro_w,real dx,real dy,real dz,int i1,int i2,int j1,int j2,int k1,int k2,real tstep, real cs2, int bcw) ;

void ic(real ***t ,real ***p ,real ***u, real ***v,real ***w ,real *ro_u,real *ro_w,real dx,real dy,real  dz, real deltaU, int i1,int i2,int j1,int j2,int k1,int k2,int bcW,int nx,int ny,int nz,real *x0,real *y0,real *z0, real *deltaTheta,real *deltaV,real thetaBar, real *rx, real *ry, real *rz,real g);

//void ic(float **t ,float **p ,float **u ,float **w ,float *ro_u,float *ro_w,float dx,float  dz,int i1,int i2,int k1,int k2,int bcW,int nx,int nz,float *x0,float *z0, float *deltaTheta,float  thetaBar,float rx,float rz,float g);

void bc4T(real ***t1, real ***t2, int i1,int i2,int j1,int j2,int k1,int k2,int bcw);
void bc(real ***u,real ***v,real ***w, int i1,int i2,int j1,int j2,int k1,int k2,int bcw);



void advection(real ***q2,real ***u,real ***v,real ***w, real ***q1, real ***uu,real ***vv, real ***ww,real dx,real dy,real dz,real dt,int i1,int i2,int j1,int j2,int k1,int k2,int nxdim,int nydim,int nzdim,real tstep, char advection_type);

void stats(float **q,int i1,int i2,int j1,int j2,int n,float *qmax, float *qmin);
void update(real ****cube1, real ****cube2);
void plot1d(float *q,int i1,int i2,int nx,int n,float qmax,int overlay,float *qtrue,char *title,char *name);
void sfc(float **q,int nx,int ny,float time,float angh,float angv,char *label,char *name);
void contr(float **z,int nx,int ny,float cint,float time,char *label,int colors,int pltzero,int nestX1,int nestX2,int nestY1,int nestY2,char *name);
//void getTakacsError(float **q1,float dx,float dy,int i1,int i2,int j1,int j2,int nx,int ny,float x0,float y0,float r);

void diffusion(real ***t,real ***u,real ***v,real ***w,real ***tt,real ***uu,real ***vv,real ***ww,real dx2,real dy2,real dz2,int i1,int i2,int j1,int j2,int k1,int k2,real dt,real tstep,real Ku,real Kv,real Kw,real Kt);


void putfield(char *name,float datatime,float *field,int nx,int ny,int nz);
// end  of prototypes

