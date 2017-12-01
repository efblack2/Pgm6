#include "real.h"
#include "constant.cuh"
#include "macros.h"

//u = u[idx3D(level,row,col,nydim,nxdim+1)] 
//v = v[idx3D(level,row,col,nydim+1,nxdim)] 
//w = w[idx3D(level,row,col,nydim,nxdim)] 
//t = t[idx3D(level,row,col,nydim,nxdim)] 
//p = p[idx3D(level-kB,row-jS,col-iW,ny,nx)] 

//#include<stdio.h>
__global__
void pgf2(real *__restrict__  p3,const real *__restrict__  u,const real *__restrict__  v,const real *__restrict__  w,const real *__restrict__  p1,const real *ro_u,const real *ro_w) 
{
    const int col = blockIdx.x*blockDim.x+threadIdx.x + iW;
    const int row = blockIdx.y*blockDim.y+threadIdx.y + jS;
    const real tempz = tstep*cs2/dz;
    
    real currentW = w[idx3D(kB,row,col,nydim,nxdim)];
    for(int level=kB; ( col<=iE  &&  row<=jN && level<=kT ); ++level) {
        real tempx = ro_u[level]*tstep*cs2/dx;
        real tempy = ro_u[level]*tstep*cs2/dy;
        real topW = w[idx3D(level+1,row,col,nydim,nxdim)];
        p3[idx3D(level-kB,row-jS,col-iW,ny,nx)] = p1[idx3D(level-kB,row-jS,col-iW,ny,nx)] - (
                                                  tempx * ( u[idx3D(level,row,col+1,nydim,nxdim+1)] - u[idx3D(level,row,col,nydim,nxdim+1)] ) +  
                                                  tempy * ( v[idx3D(level,row+1,col,nydim+1,nxdim)] - v[idx3D(level,row,col,nydim+1,nxdim)] ) +  
                                                  tempz * ( ro_w[level+1]*topW - ro_w[level]*currentW ) 
                                                  );
        currentW=topW;                                              
    } // end for //    
} // end of pgf2() //
