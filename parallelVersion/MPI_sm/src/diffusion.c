
#include "real.h"
void diffusion(real ***t,real ***u,real ***v,real ***w,real ***restrict tt,real ***restrict uu,real ***restrict vv,real ***restrict ww,real dx2,real dy2,real dz2,int i1,int i2,int j1,int j2,int k1,int k2,real dt,real tstep,real Ku,real Kv,real Kw,real Kt) 
{ 

    //#pragma omp for nowait
    for (int level=k1; level<=k2; ++level) {
        for (int row=j1;row<=j2; ++row) {
            for (int col=i1;col<=i2; ++col) {
                u[level][row][col] += tstep*Ku*( 
                                    (uu[level][row][col+1] -2*uu[level][row][col] + uu[level][row][col-1]) / (dx2) + 
                                    (uu[level][row+1][col] -2*uu[level][row][col] + uu[level][row-1][col]) / (dy2) +
                                    (uu[level+1][row][col] -2*uu[level][row][col] + uu[level-1][row][col]) / (dz2) );
            } // end for //
        } // end for //
    }  // end for //
    
    //#pragma omp for nowait
    for (int level=k1; level<=k2; ++level) {
        for (int row=j1;row<=j2; ++row) {
            for (int col=i1;col<=i2; ++col) {
                v[level][row][col] += tstep*Kv*( 
                                    (vv[level][row][col+1] -2*vv[level][row][col] + vv[level][row][col-1]) / (dx2) + 
                                    (vv[level][row+1][col] -2*vv[level][row][col] + vv[level][row-1][col]) / (dy2) +
                                    (vv[level+1][row][col] -2*vv[level][row][col] + vv[level-1][row][col]) / (dz2) );
            } // end for //
        } // end for //
    }  // end for //

    //#pragma omp for nowait
    for (int level=k1; level<=k2; ++level) {
        for (int row=j1;row<=j2; ++row) {
            for (int col=i1;col<=i2; ++col) {
                w[level][row][col] += tstep*Kw*( 
                                    (ww[level][row][col+1] -2*ww[level][row][col] + ww[level][row][col-1]) / (dx2) + 
                                    (ww[level][row+1][col] -2*ww[level][row][col] + ww[level][row-1][col]) / (dy2) +
                                    (ww[level+1][row][col] -2*ww[level][row][col] + ww[level-1][row][col]) / (dz2) );                                    
            } // end for //
        } // end for //
    }  // end for //


    //#pragma omp for // using the implicit barrier to syncromize all thread before returning
    for (int level=k1; level<=k2; ++level) {
        for (int row=j1;row<=j2; ++row) {
            for (int col=i1;col<=i2; ++col) {
                t[level][row][col] += dt*Kt*( 
                                    (tt[level][row][col+1] -2*tt[level][row][col] + tt[level][row][col-1]) / (dx2) + 
                                    (tt[level][row+1][col] -2*tt[level][row][col] + tt[level][row-1][col]) / (dy2) +
                                    (tt[level+1][row][col] -2*tt[level][row][col] + tt[level-1][row][col]) / (dz2) );
            } // end for //
        } // end for //
    }  // end for //
} // end of diffussion() //
