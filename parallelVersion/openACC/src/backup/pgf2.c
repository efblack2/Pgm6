#include "real.h"

void pgf2(real ***restrict p3,real ***restrict u,real ***restrict v,real ***restrict w,real ***restrict p1,real *ro_u,real *ro_w,real dx,real dy,real dz,int i1,int i2,int j1,int j2,int k1,int k2,real tstep, real cs2) 
{
    #pragma acc parallel  present(p3,u,v,w,p1,ro_u,ro_w) async vector_length(128)
    {
        real temp1 = -tstep*cs2/dz;
        #pragma acc loop gang
        for(int level=k1; level<=k2; ++level) {
            real tempx = -tstep*cs2*ro_u[level]/dx;
            real tempy = -tstep*cs2*ro_u[level]/dy;
            #pragma acc loop seq
            for(int row=j1; row<=j2; ++row) {
                #pragma acc loop vector
                for(int col=i1; col<=i2; ++col) {
                    p3[level-k1][row-j1][col-i1] = p1[level-k1][row-j1][col-i1] + 
                                                ( tempx * (u[level][row][col+1] - u[level][row][col]) +  
                                                  tempy * (v[level][row+1][col] - v[level][row][col]) +  
                                                  temp1 * ( ro_w[level+1]*w[level+1][row][col] -  ro_w[level]*w[level][row][col] ));
                } // end for //
            } // end for //
        } // end for //     
    } // end of parallel region
} // end of pgf2() //
