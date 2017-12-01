#include "real.h"

void pgf1(real ***restrict u,real ***restrict v,real ***restrict w,real ***restrict t,real ***restrict p1,real *ro_u,real *ro_w,real dx,real dy,real dz,int i1,int i2,int j1,int j2,int k1,int k2,real tstep,real g,real thetaBar) 
{
    #pragma acc parallel  present(u,v,w,t,p1,ro_u,ro_w) async vector_length(128)
    {
        #pragma acc loop gang
        for(int level=k1; level<=k2; ++level) {
            real tempx = -tstep/(ro_u[level]*dx);
            #pragma acc loop seq
            for(int row=j1; row<=j2; ++row) {
                #pragma acc loop vector
                for(int col=i1+1; col<=i2; ++col) {
                    u[level][row][col] += tempx * (p1[level-k1][row-j1][col-i1]-p1[level-k1][row-j1][col-i1-1]);
                } // end for //
            } // end for //
        } // end for //
        
        
        #pragma acc loop gang
        for(int level=k1; level<=k2; ++level) {
            real tempy = -tstep/(ro_u[level]*dy);
            #pragma acc loop seq
            for(int row=j1+1; row<=j2; ++row) {
                #pragma acc loop vector
                for(int col=i1; col<=i2; ++col) {
                    v[level][row][col] += tempy * (p1[level-k1][row-j1][col-i1]-p1[level-k1][row-j1-1][col-i1]);
                } // end for //
            } // end for //
        } // end for //

        
        real temp1 = tstep*g/(2.0*thetaBar);
        #pragma acc loop gang
        for(int level=k1+1; level<=k2; ++level) {
            real tempz = -tstep/(ro_w[level]*dz);
            #pragma acc loop seq
            for(int row=j1; row<=j2; ++row) {
                #pragma acc loop vector
                for(int col=i1; col<=i2; ++col) {
                    w[level][row][col] += (tempz * (p1[level-k1][row-j1][col-i1]-p1[level-k1-1][row-j1][col-i1]) + temp1*(t[level][row][col] + t[level-1][row][col] - 2.0*thetaBar) ) ;
                } // end for //
            } // end for //
        } // end for //

    } // end of parallel region
    
} // end of pgf1() //
