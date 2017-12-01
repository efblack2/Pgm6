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
void bcEW(real *u,real *v,real *w, real*t1, real *t2)
{
    const int row   = blockIdx.x*blockDim.x+threadIdx.x + jS;
    const int level = blockIdx.y*blockDim.y+threadIdx.y + kB;
    

    //  East and West faces for u  //
    for (int i=0; i<=bcw && row<=jN && level<=kT ; ++i) { 
        u[idx3D(level,row,iW-i,nydim,nxdim+1)]   = -u[idx3D(level,row,iW+1+i,nydim,nxdim+1)];
        u[idx3D(level,row,iE+1+i,nydim,nxdim+1)] = -u[idx3D(level,row,iE-i,nydim,nxdim+1)];
    } // end for //


    //  East and West faces for v  //
    for (int i=1; i<=bcw  && row<=jN+1 && level<=kT; ++i) {        
        v[idx3D(level,row,iW-i,nydim+1,nxdim)] = v[idx3D(level,row,iW+i,nydim+1,nxdim)];
        v[idx3D(level,row,iE+i,nydim+1,nxdim)] = v[idx3D(level,row,iE-i,nydim+1,nxdim)];
    } // end for //           


    //  East and West faces for w  //
    for (int i=1; i<=bcw && row<=jN && level<=kT+1; ++i) {        
        w[idx3D(level,row,iW-i,nydim,nxdim)] = w[idx3D(level,row,iW+i,nydim,nxdim)];
        w[idx3D(level,row,iE+i,nydim,nxdim)] = w[idx3D(level,row,iE-i,nydim,nxdim)];
    } // end for //
    
    //  East and West faces for t  //
    for (int i=1; i<=bcw && row<=jN && level<=kT ; ++i) {        
        t1[idx3D(level,row,iW-i,nydim,nxdim)] = t1[idx3D(level,row,iW+i,nydim,nxdim)];            
        t1[idx3D(level,row,iE+i,nydim,nxdim)] = t1[idx3D(level,row,iE-i,nydim,nxdim)];            
        t2[idx3D(level,row,iW-i,nydim,nxdim)] = t1[idx3D(level,row,iW+i,nydim,nxdim)];
        t2[idx3D(level,row,iE+i,nydim,nxdim)] = t1[idx3D(level,row,iE-i,nydim,nxdim)];
    } // end for //

} // end of bcEW() //
