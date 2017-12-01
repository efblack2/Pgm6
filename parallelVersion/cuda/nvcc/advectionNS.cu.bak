#include "real.h"
#include "constant.cuh"
#include "macros.h"

/*
 * ======================= advection NS ====================
 * Integrate forward (advection only) by one time step.
 * ATMS 502 / CSE 566, Spring 2016
 *
 * Arguments:
 *
 *	q1	real array	values at current step
 *	q2	real array	values at next step
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
void advectionNS(real *__restrict__ tTemp, const real *__restrict__ t, const real *v) // advection in y (N-S) //
{
    const int col = blockIdx.x*blockDim.x+threadIdx.x + iW;
    const int row = blockIdx.y*blockDim.y+threadIdx.y + jS;
    real tR,tRm1;
    real vR,vRp1;

    if (advection_type == 'p') {
        __shared__ real flux[TILESIZEX][TILESIZEY+1];    
        //if(col<=iE && level<=kT) {
        for (int level=kB; row<=jN && col<=iE && level<=kT ; ++level) {	            
            tR   = t[idx3D(level,row,  col,nydim,nxdim)];
            tRm1 = t[idx3D(level,row-1,col,nydim,nxdim)];

            vR   = v[idx3D(level,row,  col,nydim+1,nxdim)];
            vRp1 = v[idx3D(level,row+1,col,nydim+1,nxdim)];
            
            real r = fabs(vR*dt/dy);
            if (vR >=0.0 ) {
                flux[threadIdx.x][threadIdx.y] = r*(tRm1 + 0.25*(1-r)*(  tR  - t[idx3D(level,row-2,col,nydim,nxdim)]) );
            } else {
                flux[threadIdx.x][threadIdx.y] = r*(-tR + 0.25*(1-r)*( t[idx3D(level,row+1,col,nydim,nxdim)] - tRm1 ));	    
            } // end if 
            
            if (threadIdx.y == blockDim.y-1  || row==jN) { 
                r = fabs(vRp1*dt/dy);
                if (vRp1 >=0.0) {
                    flux[threadIdx.x][threadIdx.y+1] = r*(tR + 0.25*(1-r)*(  t[idx3D(level,row+1,col,nydim,nxdim)]-tRm1) );
                } else {
                    flux[threadIdx.x][threadIdx.y+1] = r*(-t[idx3D(level,row+1,col,nydim,nxdim)] + 0.25*(1-r)*( t[idx3D(level,row+2,col,nydim,nxdim)] - tR ));	    
                } // end if 
            } // end if //
                      
            __syncthreads();                
            tTemp[idx3D(level,row,col,nydim,nxdim)] = tR - (flux[threadIdx.x][threadIdx.y+1]-flux[threadIdx.x][threadIdx.y] ) + dt/dy*tR*(vRp1-vR);                               
            __syncthreads();
                            
        } // end for //    
        //} // end if //
    } // end if //       

} // end of advectionNS() //
