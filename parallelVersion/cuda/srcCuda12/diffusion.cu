#include "real.h"
#include "constant.cuh"
#include "macros.h"


//u = u[idx3D(level,row,col,nydim,nxdim+1)] 
//v = v[idx3D(level,row,col,nydim+1,nxdim)] 
//w = w[idx3D(level,row,col,nydim,nxdim)] 
//t = t[idx3D(level,row,col,nydim,nxdim)] 
//p = p1[idx3D(level-kB,row-jS,col-iW,ny,nx)] 
//#include <stdio.h>
__global__
void diffusion(real *t,real *u,real *v,real *w,const real *__restrict__  tt,const real *__restrict__  uu,const real *__restrict__  vv,const real *__restrict__  ww)
{
    const int col = blockIdx.x*blockDim.x+threadIdx.x +iW;
    const int row = blockIdx.y*blockDim.y+threadIdx.y +jS;
    const real dx2=dx*dx;
    const real dy2=dy*dy;
    const real dz2=dz*dz;
    constexpr real TWO = (real) 2;
    //real bot,cur,top;

    __shared__ real planeYX[TILESIZEY][TILESIZEX];

    if(col<=iE && row<=jN ) {
        real bot=uu[idx3D(kB-1,row,col,nydim,nxdim+1)];
        real cur=uu[idx3D(kB  ,row,col,nydim,nxdim+1)];
        
        for (int level=kB;  level<=kT; ++level) {
            real top=uu[idx3D(level+1,row,col,nydim,nxdim+1)];
            planeYX[threadIdx.y][threadIdx.x] = cur;
            __syncthreads();
            
            u[idx3D(level,row,col,nydim,nxdim+1)]  += tstep*KU*(    
                                                      (( (threadIdx.x < blockDim.x-1 && col < iE )  ? planeYX[threadIdx.y][threadIdx.x+1] : uu[idx3D(level,row,col+1,nydim,nxdim+1)] ) 
                                                                                                               - TWO*cur + 
                                                       (            threadIdx.x > 0                 ? planeYX[threadIdx.y][threadIdx.x-1] : uu[idx3D(level,row,col-1,nydim,nxdim+1)] ) ) / dx2 +
                                                      (( (threadIdx.y < blockDim.y-1  && row < jN)  ? planeYX[threadIdx.y+1][threadIdx.x] : uu[idx3D(level,row+1,col,nydim,nxdim+1)] ) 
                                                                                                               - TWO*cur +                                                     
                                                       (             threadIdx.y > 0                ? planeYX[threadIdx.y-1][threadIdx.x] : uu[idx3D(level,row-1,col,nydim,nxdim+1)] ) ) / dy2 +
                                                      (top - TWO*cur + bot ) / (dz2) );
                                                       
            __syncthreads();
            bot=cur;
            cur=top;
        }  // end for //
    } // end if //        

 
    if(col<=iE && row<=jN ) {
        real bot=vv[idx3D(kB-1,row,col,nydim+1,nxdim)];
        real cur=vv[idx3D(kB  ,row,col,nydim+1,nxdim)];
        
        for (int level=kB;  level<=kT; ++level) {
            real top=vv[idx3D(level+1,row,col,nydim+1,nxdim)];
            planeYX[threadIdx.y][threadIdx.x] = cur;
             __syncthreads();
            
            v[idx3D(level,row,col,nydim+1,nxdim)]  += tstep*KV*(    
                                                      (( (threadIdx.x < blockDim.x-1 && col < iE )  ? planeYX[threadIdx.y][threadIdx.x+1] : vv[idx3D(level,row,col+1,nydim+1,nxdim)] ) 
                                                                                                               - TWO*cur + 
                                                       (            threadIdx.x > 0                 ? planeYX[threadIdx.y][threadIdx.x-1] : vv[idx3D(level,row,col-1,nydim+1,nxdim)] ) ) / dx2 +
                                                      (( (threadIdx.y < blockDim.y-1  && row < jN)  ? planeYX[threadIdx.y+1][threadIdx.x] : vv[idx3D(level,row+1,col,nydim+1,nxdim)] ) 
                                                                                                               - TWO*cur +                                                     
                                                       (             threadIdx.y > 0                ? planeYX[threadIdx.y-1][threadIdx.x] : vv[idx3D(level,row-1,col,nydim+1,nxdim)] ) ) / dy2 +
                                                      (top - TWO*cur + bot ) / (dz2) );
                                                       
            __syncthreads();
            bot=cur;
            cur=top;
        }  // end for //
    } // end if //        



    if(col<=iE && row<=jN ) {
        real bot=ww[idx3D(kB-1,row,col,nydim,nxdim)];
        real cur=ww[idx3D(kB  ,row,col,nydim,nxdim)];
        
        for (int level=kB;  level<=kT; ++level) {
            real top=ww[idx3D(level+1,row,col,nydim,nxdim)];
            planeYX[threadIdx.y][threadIdx.x] = cur;
             __syncthreads();
            
            w[idx3D(level,row,col,nydim,nxdim)]  += tstep*KW*(    
                                                      (( (threadIdx.x < blockDim.x-1 && col < iE )  ? planeYX[threadIdx.y][threadIdx.x+1] : ww[idx3D(level,row,col+1,nydim,nxdim)] ) 
                                                                                                               - TWO*cur + 
                                                       (            threadIdx.x > 0                 ? planeYX[threadIdx.y][threadIdx.x-1] : ww[idx3D(level,row,col-1,nydim,nxdim)] ) ) / dx2 +
                                                      (( (threadIdx.y < blockDim.y-1  && row < jN)  ? planeYX[threadIdx.y+1][threadIdx.x] : ww[idx3D(level,row+1,col,nydim,nxdim)] ) 
                                                                                                               - TWO*cur +                                                     
                                                       (             threadIdx.y > 0                ? planeYX[threadIdx.y-1][threadIdx.x] : ww[idx3D(level,row-1,col,nydim,nxdim)] ) ) / dy2 +
                                                      (top - TWO*cur + bot ) / (dz2) );
                                                       
            __syncthreads();
            bot=cur;
            cur=top;
        }  // end for //
    } // end if //        

    if(col<=iE && row<=jN ) {    
        real bot=tt[idx3D(kB-1,row,col,nydim,nxdim)];
        real cur=tt[idx3D(kB  ,row,col,nydim,nxdim)];
        
        for (int level=kB; level<=kT; ++level) {
            real top=tt[idx3D(level+1,row,col,nydim,nxdim)];
            planeYX[threadIdx.y][threadIdx.x] = cur;
             __syncthreads();
            
            t[idx3D(level,row,col,nydim,nxdim)]  += dt*KTemp*(    
                                                      (( (threadIdx.x < blockDim.x-1 && col < iE )  ? planeYX[threadIdx.y][threadIdx.x+1] : tt[idx3D(level,row,col+1,nydim,nxdim)] ) 
                                                                                                               - TWO*cur + 
                                                       (            threadIdx.x > 0                 ? planeYX[threadIdx.y][threadIdx.x-1] : tt[idx3D(level,row,col-1,nydim,nxdim)] ) ) / dx2 +
                                                      (( (threadIdx.y < blockDim.y-1  && row < jN)  ? planeYX[threadIdx.y+1][threadIdx.x] : tt[idx3D(level,row+1,col,nydim,nxdim)] ) 
                                                                                                               - TWO*cur +                                                     
                                                       (             threadIdx.y > 0                ? planeYX[threadIdx.y-1][threadIdx.x] : tt[idx3D(level,row-1,col,nydim,nxdim)] ) ) / dy2 +
                                                      (top - TWO*cur + bot ) / (dz2) );                                                       
            __syncthreads();
            bot=cur;
            cur=top;
        }  // end for //
    } // end if //     


/*

    if(col<=iE && row<=jN ) {
        real bot=uu[idx3D(kB-1,row,col,nydim,nxdim+1)];
        real cur=uu[idx3D(kB,row,col,nydim,nxdim+1)];
        for (int level=kB; level<=kT; ++level) {
            real top=uu[idx3D(level+1,row,col,nydim,nxdim+1)];
            u[idx3D(level,row,col,nydim,nxdim+1)]  += tstep*KU*( 
                                                        (uu[idx3D(level,row,col+1,nydim,nxdim+1)] - TWO*cur + uu[idx3D(level,row,col-1,nydim,nxdim+1)] ) / (dx2) + 
                                                        (uu[idx3D(level,row+1,col,nydim,nxdim+1)] - TWO*cur + uu[idx3D(level,row-1,col,nydim,nxdim+1)] ) / (dy2) +
                                                        (top - TWO*cur + bot ) / (dz2) );
            bot=cur;
            cur=top;
        }  // end for //
    } // end if //    


    if(col<=iE && row<=jN ) {    
        real bot=vv[idx3D(kB-1,row,col,nydim+1,nxdim)];    
        real cur=vv[idx3D(kB,row,col,nydim+1,nxdim)];    
        for (int level=kB; level<=kT; ++level) {
            real top=vv[idx3D(level+1,row,col,nydim+1,nxdim)];    
            v[idx3D(level,row,col,nydim+1,nxdim)] += tstep*KV*( 
                                                        (vv[idx3D(level,row,col+1,nydim+1,nxdim)] - TWO*cur + vv[idx3D(level,row,col-1,nydim+1,nxdim)]) / dx2 + 
                                                        (vv[idx3D(level,row+1,col,nydim+1,nxdim)] - TWO*cur + vv[idx3D(level,row-1,col,nydim+1,nxdim)]) / dy2 +
                                                        (top - TWO*cur + bot ) / dz2 );
            bot=cur;
            cur=top;
        }  // end for //
    } // end if //    


    if(col<=iE && row<=jN ) {    
        real bot=ww[idx3D(kB-1,row,col,nydim,nxdim)];
        real cur=ww[idx3D(kB,row,col,nydim,nxdim)];
        for (int level=kB; level<=kT; ++level) {
            real top=ww[idx3D(level+1,row,col,nydim,nxdim)];
            w[idx3D(level,row,col,nydim,nxdim)] += tstep*KW*( 
                                                        (ww[idx3D(level,row,col+1,nydim,nxdim)]  - TWO*cur + ww[idx3D(level,row,col-1,nydim,nxdim)]) / dx2 + 
                                                        (ww[idx3D(level,row+1,col,nydim,nxdim)]  - TWO*cur + ww[idx3D(level,row-1,col,nydim,nxdim)]) / dy2 +
                                                        (top - TWO*cur + bot ) / dz2 );
            bot=cur;
            cur=top;
        }  // end for //
    } // end if //    




    if(col<=iE && row<=jN ) {    
        real bot=tt[idx3D(kB-1,row,col,nydim,nxdim)];
        real cur=tt[idx3D(kB,row,col,nydim,nxdim)];
        for (int level=kB; level<=kT; ++level) {
            real top=tt[idx3D(level+1,row,col,nydim,nxdim)];
            t[idx3D(level,row,col,nydim,nxdim)] += dt*KTemp*( 
                                                        (tt[idx3D(level,row,col+1,nydim,nxdim)]  - TWO*cur + tt[idx3D(level,row,col-1,nydim,nxdim)]) / dx2 + 
                                                        (tt[idx3D(level,row+1,col,nydim,nxdim)]  - TWO*cur + tt[idx3D(level,row-1,col,nydim,nxdim)]) / dy2 +
                                                        (top - TWO*cur + bot ) / dz2 );
            bot=cur;
            cur=top;
        }  // end for //
    } // end if //    


*/    
} // end of diffussion() //
