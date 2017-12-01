#include <math.h>
#include "real.h"

real Density(real z,real g,real thetaBar)
{
    const real  cp = 1004.;
    const real  Rd = 287.038;
    const real  P0 = 1.0e5;
    float T,P;
    T=thetaBar - g*z/cp;
    P=P0*pow((T/thetaBar),(cp/Rd));
    //printf("\t\tz, T, P: %10.4f %10.4f %10.4f %10.4f\n",z,T,P,P/(Rd*T)); 
    return P/(Rd*T);
} // end of Density() //

void getDensity(real *ro_u,real *ro_w,int kB,int kT,real dz,real g,real thetaBar)
{
    real  z=0.0;
    for (int row=kB; row<=kT; ++row){
        z = dz*(row-kB);
        ro_w[row] = Density(z,g,thetaBar);
        z+=0.5*dz;
        ro_u[row] = Density(z,g,thetaBar);
    } // end for //
    /*
    for (int row=kB+1; row<=kT; ++row){
        ro_w[row] = 0.5*(ro_u[row]+ro_u[row-1]);
    } // end for //
    */
    
    //z+=0.5*dz;
    ro_w[kT+1] = 0.0; ;//Density(z,g,thetaBar);
} // end of getDensity() //
