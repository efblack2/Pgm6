#include "real.h"
#include "constant.cuh"
#include "macros.h"

/*
 * ======================= advection ====================
 * Integrate forward (advection only) by one time step.
 * ATMS 502 / CSE 566, Spring 2016
 *
 * Arguments:
 *
 *	q1	real array	values at current step
 *	q2	real array	values at next step
 *	c	real		true speed of wave
 *	dx	real		grid spacing
 *	iW,iE	integers	indices bounding array data
 *	nx	integer		number of grid points
 *	advection_type
 *              char 		if 'L', linear advection;
 *				otherwise, nonlinear
 */
//u = u[idx3D(level,row,col,nydim,nxdim+1)] 
//v = v[idx3D(level,row,col,nydim+1,nxdim)] 
//w = w[idx3D(level,row,col,nydim,nxdim)] 
#include <stdio.h>
__global__
void advecUVW(real *u,real *v,real *w, const real *__restrict__  uu,const real *__restrict__  vv,const real *__restrict__  ww)

{

    const int col = blockIdx.x*blockDim.x+threadIdx.x + iW;
    const int row = blockIdx.y*blockDim.y+threadIdx.y + jS;
    
    __shared__ real planeYX[TILESIZEY][TILESIZEX];

    constexpr real FOUR = (real) 4;
    const real tempx = tstep/(FOUR*dx);
    const real tempy = tstep/(FOUR*dy);
    const real tempz = tstep/(FOUR*dz);
    real bot,cur;  


//    if () {
    bot=uu[idx3D(kB-1,row,col,nydim,nxdim+1)];
    cur=uu[idx3D(kB  ,row,col,nydim,nxdim+1)];
    for (int level=kB; col<=iE+1 && row<=jN && level<=kT; ++level) {	
        real top=uu[idx3D(level+1,row,col,nydim,nxdim+1)];
        planeYX[threadIdx.y][threadIdx.x] = cur;
        __syncthreads();            
        u[idx3D(level,row,col,nydim,nxdim+1)] -= 
        tempx *(
              (cur+(             threadIdx.x > 0            ? planeYX[threadIdx.y][threadIdx.x-1] : uu[idx3D(level,row,col-1,nydim,nxdim+1)]))*
              (cur-(             threadIdx.x > 0            ? planeYX[threadIdx.y][threadIdx.x-1] : uu[idx3D(level,row,col-1,nydim,nxdim+1)])) +
              (( (threadIdx.x < blockDim.x-1 && col <= iE ) ? planeYX[threadIdx.y][threadIdx.x+1] : uu[idx3D(level,row,col+1,nydim,nxdim+1)])+cur)*  
              (( (threadIdx.x < blockDim.x-1 && col <= iE ) ? planeYX[threadIdx.y][threadIdx.x+1] : uu[idx3D(level,row,col+1,nydim,nxdim+1)]) - cur)
               ) +
        tempy *(
              (vv[idx3D(level,row  ,col,nydim+1,nxdim)] + vv[idx3D(level,row  ,col-1,nydim+1,nxdim)]) * 
              (cur-(             threadIdx.y > 0            ? planeYX[threadIdx.y-1][threadIdx.x] : uu[idx3D(level,row-1,col,nydim,nxdim+1)]) ) +
              (vv[idx3D(level,row+1,col,nydim+1,nxdim)] + vv[idx3D(level,row+1,col-1,nydim+1,nxdim)]) *
              (( (threadIdx.y < blockDim.y-1  && row < jN)  ? planeYX[threadIdx.y+1][threadIdx.x] : uu[idx3D(level,row+1,col,nydim,nxdim+1)]) - cur)
              ) +                    
        tempz *(
              (ww[idx3D(level  ,row,col,nydim,nxdim)] + ww[idx3D(level  ,row,col-1,nydim,nxdim)])* 
              (cur - bot)  +
              (ww[idx3D(level+1,row,col,nydim,nxdim)] + ww[idx3D(level+1,row,col-1,nydim,nxdim)])* 
              (top - cur)
              );
        __syncthreads();
        bot=cur;
        cur=top;
    } // end for //
//    } // end if //


//    if () {
    bot=vv[idx3D(kB-1,row,col,nydim+1,nxdim)];
    cur=vv[idx3D(kB  ,row,col,nydim+1,nxdim)];
    for (int level=kB; col<=iE && row<=jN+1 && level<=kT; ++level) {
        real top = vv[idx3D(level+1,row,col,nydim+1,nxdim)];
        planeYX[threadIdx.y][threadIdx.x] = cur;
        __syncthreads();
        v[idx3D(level,row,col,nydim+1,nxdim)] -= 
        tempx *(
            (uu[idx3D(level,row,col  ,nydim,nxdim+1)] + uu[idx3D(level,row-1,col  ,nydim,nxdim+1)]) * 
            (cur-(            threadIdx.x > 0             ? planeYX[threadIdx.y][threadIdx.x-1] : vv[idx3D(level,row,col-1,nydim+1,nxdim)]) ) + 
            (uu[idx3D(level,row,col+1,nydim,nxdim+1)] + uu[idx3D(level,row-1,col+1,nydim,nxdim+1)]) * 
            (( (threadIdx.x < blockDim.x-1 && col < iE )  ? planeYX[threadIdx.y][threadIdx.x+1] : vv[idx3D(level,row,col+1,nydim+1,nxdim)]) -cur)
            )+
        tempy *(
            (cur+(             threadIdx.y > 0              ? planeYX[threadIdx.y-1][threadIdx.x] : vv[idx3D(level,row-1,col,nydim+1,nxdim)])  )* 
            (cur-(             threadIdx.y > 0              ? planeYX[threadIdx.y-1][threadIdx.x] : vv[idx3D(level,row-1,col,nydim+1,nxdim)]) ) + 
            ( ( (threadIdx.y < blockDim.y-1  && row <= jN)  ? planeYX[threadIdx.y+1][threadIdx.x] : vv[idx3D(level,row+1,col,nydim+1,nxdim)])  + cur) * 
            ( ( (threadIdx.y < blockDim.y-1  && row <= jN)  ? planeYX[threadIdx.y+1][threadIdx.x] : vv[idx3D(level,row+1,col,nydim+1,nxdim)])  - cur)
            )+
        tempz *(
            (ww[idx3D(level  ,row,col,nydim,nxdim)] + ww[idx3D(level  ,row-1,col,nydim,nxdim)]) * 
            (cur - bot) + 
            (ww[idx3D(level+1,row,col,nydim,nxdim)] + ww[idx3D(level+1,row-1,col,nydim,nxdim)]) * 
            (top - cur)
            );
        __syncthreads();
        bot=cur;
        cur=top;
    } // end for //
//    } // end if //


//    if ( ) {
    bot=ww[idx3D(kB-1,row,col,nydim,nxdim)];
    cur=ww[idx3D(kB  ,row,col,nydim,nxdim)];
    for (int level=kB; col<=iE && row<=jN && level<=kT+1; ++level) {	
        real top = ww[idx3D(level+1,row,col,nydim,nxdim)];
        planeYX[threadIdx.y][threadIdx.x] = cur;
        __syncthreads();           
        w[idx3D(level,row,col,nydim,nxdim)] -= 
        tempx *(
            (uu[idx3D(level,row,col  ,nydim,nxdim+1)] + uu[idx3D(level-1,row,col  ,nydim,nxdim+1)])  * 
            (cur - (threadIdx.x > 0                     ? planeYX[threadIdx.y][threadIdx.x-1]:ww[idx3D(level,row,col-1,nydim,nxdim)])  ) +
            (uu[idx3D(level,row,col+1,nydim,nxdim+1)] + uu[idx3D(level-1,row,col+1,nydim,nxdim+1)])  * 
            (((threadIdx.x < blockDim.x-1 && col < iE ) ? planeYX[threadIdx.y][threadIdx.x+1]:ww[idx3D(level,row,col+1,nydim,nxdim)]) - cur)
               )+
        tempy *(
            (vv[idx3D(level,row  ,col,nydim+1,nxdim)] + vv[idx3D(level-1,row  ,col,nydim+1,nxdim)])  * 
            (cur - (threadIdx.y > 0                       ? planeYX[threadIdx.y-1][threadIdx.x]:ww[idx3D(level,row-1,col,nydim,nxdim)]) ) + 
            (vv[idx3D(level,row+1,col,nydim+1,nxdim)] + vv[idx3D(level-1,row+1,col,nydim+1,nxdim)])  * 
            (( (threadIdx.y < blockDim.y-1  && row < jN)  ? planeYX[threadIdx.y+1][threadIdx.x]:ww[idx3D(level,row+1,col,nydim,nxdim)]) - cur)
               ) +
        tempz *(
            (cur + bot ) * 
            (cur - bot)  +
            (top + cur)  * 
            (top - cur)
               );
        __syncthreads();
        bot=cur;
        cur=top;
    } // end for //
//    } // end if //
} // end of advecUVW() //    
