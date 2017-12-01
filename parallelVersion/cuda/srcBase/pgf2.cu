#include "real.h"

void pgf2(real ***__restrict__  p3,real ***__restrict__  u,real ***__restrict__  v,real ***__restrict__  w,real ***__restrict__  t,real ***__restrict__  p1,real *ro_u,real *ro_w,real dx,real dy,real dz,int iW,int iE,int jS,int jN,int kB,int kT,real tstep,real cs2, int bcw) 
{

    real tempz = -tstep*cs2/dz;
    for(int level=kB; level<=kT; ++level) {
        real tempx = -tstep*cs2*ro_u[level]/dx;
        real tempy = -tstep*cs2*ro_u[level]/dy;
        for(int row=jS; row<=jN; ++row) {
            for(int col=iW; col<=iE; ++col) {
                p3[level-kB][row-jS][col-iW] = p1[level-kB][row-jS][col-iW] + 
                                            ( tempx * (u[level][row][col+1] - u[level][row][col]) +  
                                              tempy * (v[level][row+1][col] - v[level][row][col]) +  
                                              tempz * ( ro_w[level+1]*w[level+1][row][col] -  ro_w[level]*w[level][row][col] ));
            } // end for //
        } // end for //
    } // end for //    
    
} // end of pgf2() //
