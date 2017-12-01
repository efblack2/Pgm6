#include "real.h"
void diffusion(real ***t,real ***u,real ***v,real ***w,real ***__restrict__  tt,real ***__restrict__  uu,real ***__restrict__  vv,real ***__restrict__  ww,real dx2,real dy2,real dz2,int iW,int iE,int jS,int jN,int kB,int kT,real dt,real tstep,real Ku,real Kv,real Kw,real Kt) 
{
   
    for (int level=kB; level<=kT; ++level) {
        for (int row=jS;row<=jN;++row) {
            for (int col=iW;col<=iE;++col) {
                u[level][row][col] += tstep*Ku*( 
                                    (uu[level][row][col+1] -2*uu[level][row][col] + uu[level][row][col-1]) / (dx2) + 
                                    (uu[level][row+1][col] -2*uu[level][row][col] + uu[level][row-1][col]) / (dy2) +
                                    (uu[level+1][row][col] -2*uu[level][row][col] + uu[level-1][row][col]) / (dz2) );
            } // end for //
        } // end for //
    }  // end for //

    for (int level=kB; level<=kT; ++level) {
        for (int row=jS;row<=jN;++row) {
            for (int col=iW;col<=iE;++col) {
                v[level][row][col] += tstep*Kv*( 
                                    (vv[level][row][col+1] -2*vv[level][row][col] + vv[level][row][col-1]) / (dx2) + 
                                    (vv[level][row+1][col] -2*vv[level][row][col] + vv[level][row-1][col]) / (dy2) +
                                    (vv[level+1][row][col] -2*vv[level][row][col] + vv[level-1][row][col]) / (dz2) );
            } // end for //
        } // end for //
    }  // end for //

    for (int level=kB; level<=kT; ++level) {
        for (int row=jS;row<=jN;++row) {
            for (int col=iW;col<=iE;++col) {
                w[level][row][col] += tstep*Kw*( 
                                    (ww[level][row][col+1] -2*ww[level][row][col] + ww[level][row][col-1]) / (dx2) + 
                                    (ww[level][row+1][col] -2*ww[level][row][col] + ww[level][row-1][col]) / (dy2) +
                                    (ww[level+1][row][col] -2*ww[level][row][col] + ww[level-1][row][col]) / (dz2) );                                    
            } // end for //
        } // end for //
    }  // end for //

    for (int level=kB; level<=kT; ++level) {
        for (int row=jS;row<=jN;++row) {
            for (int col=iW;col<=iE;++col) {
                t[level][row][col] += dt*Kt*( 
                                    (tt[level][row][col+1] -2*tt[level][row][col] + tt[level][row][col-1]) / (dx2) + 
                                    (tt[level][row+1][col] -2*tt[level][row][col] + tt[level][row-1][col]) / (dy2) +
                                    (tt[level+1][row][col] -2*tt[level][row][col] + tt[level-1][row][col]) / (dz2) );
            } // end for //
        } // end for //
    }  // end for //
} // end of diffussion() //
