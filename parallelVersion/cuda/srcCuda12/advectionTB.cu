#include "real.h"
#include "constant.cuh"
#include "macros.h"

/*
 * ======================= advection TB ====================
 * Integrate forward (advection only) by one time step.
 * ATMS 502 / CSE 566, Spring 2016
 *
 * Arguments:
 *
 *	t1	real array	values at current step
 *	t2	real array	values at next step
 *	c	real		true speed of wave
 *	dx	real		grid spacing
 *	dt	real		time step
 *	iW,iE	integers	indices bounding array data
 *	nx	integer		number of grid points
 *	advection_type
 *              char 		if 'L', linear advection;
 *				otherwise, nonlinear
 */
//u = u[idx3D(level,row,col,nydim,nxdim+1)] 
//v = v[idx3D(level,row,col,nydim+1,nxdim)] 
//w = w[idx3D(level,row,col,nydim,nxdim)] 
//t = t[idx3D(level,row,col,nydim,nxdim)] 
//#include <stdio.h>
__global__
void advectionTB(real *__restrict__ tTemp, const real *__restrict__ t, const real *w)
{
    const int col   = blockIdx.x*blockDim.x+threadIdx.x + iW;
    const int level = blockIdx.y*blockDim.y+threadIdx.y + kB;
    real tR,tRm1;
    real wR,wRp1;

    if (advection_type == 'p') {
        __shared__ real flux[TILESIZEX][TILESIZEY+1];    
        //if(col<=iE && level<=kT) {
        for (int row=jS; row<=jN && col<=iE && level<=kT ; ++row) {	
            tR   = t[idx3D(level,row,col,nydim,nxdim)];
            tRm1 = t[idx3D(level-1,row,col,nydim,nxdim)];

            wR   = w[idx3D(level,row,col,nydim,nxdim)];
            wRp1 = w[idx3D(level+1,row,col,nydim,nxdim)];
            
            real r = fabs(wR*dt/dz);
            if (wR >=0.0 ) {
                flux[threadIdx.x][threadIdx.y] = r*(tRm1 + 0.25*(1-r)*(  tR  - t[idx3D(level-2,row,col,nydim,nxdim)]) );
            } else {
                flux[threadIdx.x][threadIdx.y] = r*(-tR + 0.25*(1-r)*( t[idx3D(level+1,row,col,nydim,nxdim)] - tRm1 ));	    
            } // end if 
            
            if (threadIdx.y == blockDim.y-1  || level==kT) { 
                r = fabs(wRp1*dt/dz);
                if (wRp1 >=0.0) {
                    flux[threadIdx.x][threadIdx.y+1] = r*(tR + 0.25*(1-r)*(  t[idx3D(level+1,row,col,nydim,nxdim)]  - tRm1) );
                } else {
                    flux[threadIdx.x][threadIdx.y+1] = r*(-t[idx3D(level+1,row,col,nydim,nxdim)] + 0.25*(1-r)*( t[idx3D(level+2,row,col,nydim,nxdim)] - tR ));	    
                } // end if 
            } // end if //
                      
            __syncthreads();
            tTemp[idx3D(level,row,col,nydim,nxdim)] = tR-(flux[threadIdx.x][threadIdx.y+1]-flux[threadIdx.x][threadIdx.y]) + dt/dz*tR*(wRp1-wR);
            __syncthreads();
            
        } // end for //    
        //} // end if //
    } // end if //       
} // end of advectionTB() //
