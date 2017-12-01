#include "real.h"
#include "constant.cuh"
#include "macros.h"

/*
 * ======================= advection EW ====================
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
void advectionEW(real *__restrict__ tTemp, const real *__restrict__ t, const real *u) // advection in x (E-W) //
{
    const int col = blockIdx.x*blockDim.x+threadIdx.x + iW;
    const int row = blockIdx.y*blockDim.y+threadIdx.y + jS;
    real tR,tRm1;
    real uR,uRp1;

    if (advection_type == 'p') {
        __shared__ real flux[TILESIZEY][TILESIZEX+1];    
        //if(col<=iE && level<=kT) {
        for (int level=kB; row<=jN && col<=iE && level<=kT ; ++level) {	            
            tR   = t[idx3D(level,row,col  ,nydim,nxdim)];
            tRm1 = t[idx3D(level,row,col-1,nydim,nxdim)];

            uR   = u[idx3D(level,row,col  ,nydim,nxdim+1)];
            uRp1 = u[idx3D(level,row,col+1,nydim,nxdim+1)];
            
            real r = fabs(uR*dt/dx);
            if (uR >=0.0 ) {
                flux[threadIdx.y][threadIdx.x] = r*(tRm1 + 0.25*(1-r)*(  tR  - t[idx3D(level,row,col-2,nydim,nxdim)]) );
            } else {
                flux[threadIdx.y][threadIdx.x] = r*(-tR + 0.25*(1-r)*( t[idx3D(level,row,col+1,nydim,nxdim)] - tRm1 ));	    
            } // end if 
            
            if (threadIdx.x == blockDim.x-1  || col==iE) { 
                r = fabs(uRp1*dt/dx);
                if (uRp1 >=0.0) {
                    flux[threadIdx.y][threadIdx.x+1] = r*(tR + 0.25*(1-r)*(  t[idx3D(level,row,col+1,nydim,nxdim)]-tRm1) );
                } else {
                    flux[threadIdx.y][threadIdx.x+1] = r*(-t[idx3D(level,row,col+1,nydim,nxdim)] + 0.25*(1-r)*( t[idx3D(level,row,col+2,nydim,nxdim)] - tR ));	    
                } // end if 
            } // end if //
                      
            __syncthreads();                
            tTemp[idx3D(level,row,col,nydim,nxdim)] = tR - (flux[threadIdx.y][threadIdx.x+1]-flux[threadIdx.y][threadIdx.x] ) + dt/dx*tR*(uRp1-uR);                               
            __syncthreads();
                            
        } // end for //    
        //} // end if //
    } // end if //       



/*
    real *__restrict__  fld;

    // advection in x (E-W) //
    fld = (real *) malloc(nxdim*sizeof(real));
    
    for (int level=kB; level<=kT; ++level) {	
        for (int row=jS; row<=jN; ++row) {	
            if (advection_type != 'p') {
                for (int col=iW; col<=iE; ++col) {
                    fld[col] = 0.5*(u[level][row][col] + u[level][row][col+1]);  
                } // end for //
            } else {
                for (int col=iW; col<=iE+1; ++col) {
                    fld[col] = u[level][row][col];
                } // end for //
            } // end if //
            advec1(q2[level][row],q1[level][row],fld,dx,dt,iW,iE,nxdim,advection_type);
        } // end for //
    } // end for //
    free(fld);
*/    

} // end of advectionEW() //
