#include "real.h"

void pgf1(real ***__restrict__  u,real ***__restrict__  v,real ***__restrict__  w,real ***__restrict__  t,real ***__restrict__  p1,real *ro_u,real *ro_w,real dx,real dy,real dz,int iW,int iE,int jS,int jN,int kB,int kT,real tstep,real g,real thetaBar, int bcw) 
{

    for(int level=kB; level<=kT; ++level) {
        real tempx = -tstep/(ro_u[level]*dx);
        for(int row=jS; row<=jN; ++row) {
            for(int col=iW+1; col<=iE; ++col) {
                u[level][row][col] += tempx * (p1[level-kB][row-jS][col-iW]-p1[level-kB][row-jS][col-iW-1]);
            } // end for //
        } // end for //
    } // end for //
    
    
    for(int level=kB; level<=kT; ++level) {
        real tempy = -tstep/(ro_u[level]*dy);
        for(int row=jS+1; row<=jN; ++row) {
            for(int col=iW; col<=iE; ++col) {
                v[level][row][col] += tempy * (p1[level-kB][row-jS][col-iW]-p1[level-kB][row-jS-1][col-iW]);
            } // end for //
        } // end for //
    } // end for //

    
    real temp1 = tstep*g/(2.0*thetaBar);
    for(int level=kB+1; level<=kT; ++level) {
        real tempz = -tstep/(ro_w[level]*dz);
        for(int row=jS; row<=jN; ++row) {
            for(int col=iW; col<=iE; ++col) {
                w[level][row][col] += (tempz * (p1[level-kB][row-jS][col-iW]-p1[level-kB-1][row-jS][col-iW]) + temp1*(t[level][row][col] + t[level-1][row][col] - 2.0*thetaBar) ) ;
            } // end for //
        } // end for //
    } // end for //
    
} // end of pgf1() //
