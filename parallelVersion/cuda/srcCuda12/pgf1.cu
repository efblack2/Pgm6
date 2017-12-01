#include "real.h"
#include "constant.cuh"
#include "macros.h"

//u = u[idx3D(level,row,col,nydim,nxdim+1)] 
//v = v[idx3D(level,row,col,nydim+1,nxdim)] 
//w = w[idx3D(level,row,col,nydim,nxdim)] 
//t = t[idx3D(level,row,col,nydim,nxdim)] 
//p = p1[idx3D(level-kB,row-jS,col-iW,ny,nx)] 

__global__

void pgf1(real *__restrict__  u,real *__restrict__  v,real *__restrict__  w,const real *__restrict__  t,const real *__restrict__  p1,const real *ro_u,const real *ro_w) 
{
    int col;
    int row;
    
    col = blockIdx.x*blockDim.x+threadIdx.x + iW+1;
    row = blockIdx.y*blockDim.y+threadIdx.y + jS;
    for(int level=kB; col<=iE  && row<=jN && level<=kT; ++level) {
        real tempx = tstep/(ro_u[level]*dx);
        u[idx3D(level,row,col,nydim,nxdim+1)] -= tempx * (p1[idx3D(level-kB,row-jS,col-iW,ny,nx)] - p1[idx3D(level-kB,row-jS,col-iW-1,ny,nx)] );
    } // end for //

    col = blockIdx.x*blockDim.x+threadIdx.x + iW ;
    row = blockIdx.y*blockDim.y+threadIdx.y + jS+1;
    for(int level=kB; col<=iE && row<=jN && level<=kT; ++level) {
        real tempy = tstep/(ro_u[level]*dy);
        v[idx3D(level,row,col,nydim+1,nxdim)] -= tempy * ( p1[idx3D(level-kB,row-jS,col-iW,ny,nx)] - p1[idx3D(level-kB,row-jS-1,col-iW,ny,nx)] );
    } // end for //

    col = blockIdx.x*blockDim.x+threadIdx.x + iW;
    row = blockIdx.y*blockDim.y+threadIdx.y + jS;
    
    real tempt = tstep*g/(2.0*thetaBar);
    real bottomP1 = p1[idx3D(0,row-jS,col-iW,ny,nx)];    
    real bottomT  =  t[idx3D(kB,row,col,nydim,nxdim)];  
    for(int level=kB+1; col<=iE && row<=jN && level<=kT; ++level) {
        real currentP1 = p1[idx3D(level-kB,row-jS,col-iW,ny,nx)];
        real currentT = t[idx3D(level,row,col,nydim,nxdim)] ;
        real tempz = tstep/(ro_w[level]*dz);
        w[idx3D(level,row,col,nydim,nxdim)] += (  tempt * ( currentT + bottomT ) - tstep*g - tempz * ( currentP1 - bottomP1 )      );
                                                 
        bottomP1 = currentP1;
        bottomT  = currentT;                                                 
    } // end for //
} // end of pgf1() //
