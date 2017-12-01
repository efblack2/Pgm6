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
void bcTB(real *u,real *v,real *w, real*t1, real *t2)
{

    const int col = blockIdx.x*blockDim.x+threadIdx.x + iW;
    const int row = blockIdx.y*blockDim.y+threadIdx.y + jS;
    real temp1,temp2;

    //  Top and bottom faces for u //
    if (col <= iE+1 && row <=jN) {
        temp1=u[idx3D(kB,row,col,nydim,nxdim+1)];
        temp2=u[idx3D(kT,row,col,nydim,nxdim+1)];
        for (int k=1; k<=bcw; ++k) {        
            u[idx3D(kB-k,row,col,nydim,nxdim+1)] = temp1;
            u[idx3D(kT+k,row,col,nydim,nxdim+1)] = temp2;
        } // end for //
    } // end if //

    //  Top and bottom faces for v //
    if (col <= iE && row <=jN+1) {
        temp1=v[idx3D(kB,row,col,nydim+1,nxdim)];// v[kB][row][col];
        temp2=v[idx3D(kT,row,col,nydim+1,nxdim)]; //v[kT][row][col];
        for (int k=1; k<=bcw; ++k) {        
            v[idx3D(kB-k,row,col,nydim+1,nxdim)]= temp1;
            v[idx3D(kT+k,row,col,nydim+1,nxdim)]= temp2;
        } // end for //
    } // end if //

    //  Top and bottom faces for w //
    if (col <= iE && row <=jN) {
        temp1 = t1[idx3D(kB,row,col,nydim,nxdim)];
        temp2 = t1[idx3D(kT,row,col,nydim,nxdim)];
        for (int k=0; k<=bcw; ++k) {        
            w[idx3D(kB-k,row,col,nydim,nxdim)]   = 0.0;
            w[idx3D(kT+1+k,row,col,nydim,nxdim)] = 0.0;
            t2[idx3D(kB-k,row,col,nydim,nxdim)] = t1[idx3D(kB-k,row,col,nydim,nxdim)] = temp1;            
            t2[idx3D(kT+k,row,col,nydim,nxdim)] = t1[idx3D(kT+k,row,col,nydim,nxdim)] = temp2;
        } // end for //
    } // end if //    
} // end of bcTB() //
