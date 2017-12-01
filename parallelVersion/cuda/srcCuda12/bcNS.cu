#include "real.h"
#include "constant.cuh"
#include "macros.h"

/*
 * ============================ bc =====================
 * BC sets the boundary conditions
 * ATMS 502 / CSE 566, Spring 2016
 *
 * Arguments:
 *
 *	q1	real array	values at current time step
 *	iW,iE	integers	indices bounding array data
 *	nx	integer		main array size, not including
 *				extra 'ghost' zones/points
 */
//u = u[idx3D(level,row,col,nydim,nxdim+1)] 
//v = v[idx3D(level,row,col,nydim+1,nxdim)] 
//w = w[idx3D(level,row,col,nydim,nxdim)] 
//t = t[idx3D(level,row,col,nydim,nxdim)] 

__global__
void bcNS(real *u,real *v,real *w, real*t1, real *t2)
{
    const int col   = blockIdx.x*blockDim.x+threadIdx.x + iW;
    const int level = blockIdx.y*blockDim.y+threadIdx.y + kB;
    

    //  North and South faces for u (periodic) //
    for (int j=1; j<=bcw && col<=iE+1 && level<=kT; ++j) {        
        u[idx3D(level,jS-j,col,nydim,nxdim+1)] = u[idx3D(level,jN+1-j,col,nydim,nxdim+1)];
        u[idx3D(level,jN+j,col,nydim,nxdim+1)] = u[idx3D(level,jS+j-1,col,nydim,nxdim+1)];
    } // end for //    
 

    //  North and South faces for v (periodic) //
    if ( col<=iE && level<=kT) {
        v[idx3D(level,jN+1,col,nydim+1,nxdim)] = v[idx3D(level,jS,col,nydim+1,nxdim)];
        for (int j=1; j<=bcw; ++j) {        
            v[idx3D(level,jS-j  ,col,nydim+1,nxdim)] = v[idx3D(level,jN+1-j,col,nydim+1,nxdim)];
            v[idx3D(level,jN+1+j,col,nydim+1,nxdim)] = v[idx3D(level,jS+j  ,col,nydim+1,nxdim)];
        } // end for //
    } // end if //


    //  North and South faces for w  (periodic) //
    for (int j=1; j<=bcw && col<=iE && level<=kT+1; ++j) {        
        w[idx3D(level,jS-j,col,nydim,nxdim)] = w[idx3D(level,jN+1-j,col,nydim,nxdim)];
        w[idx3D(level,jN+j,col,nydim,nxdim)] = w[idx3D(level,jS+j-1,col,nydim,nxdim)];
    } // end for //
    
    //  North and South faces for t  (periodic) //
    for (int j=1; j<=bcw && col<=iE && level<=kT; ++j) {        
        t1[idx3D(level,jS-j,col,nydim,nxdim)] = t1[idx3D(level,jN+1-j,col,nydim,nxdim)];
        t1[idx3D(level,jN+j,col,nydim,nxdim)] = t1[idx3D(level,jS+j-1,col,nydim,nxdim)];
        t2[idx3D(level,jS-j,col,nydim,nxdim)] = t1[idx3D(level,jN+1-j,col,nydim,nxdim)];
        t2[idx3D(level,jN+j,col,nydim,nxdim)] = t1[idx3D(level,jS+j-1,col,nydim,nxdim)];
    } // end for //
        
} // end of bcNS() //
