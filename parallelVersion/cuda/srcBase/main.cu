/*
 *  ATMS 502 / CSE 566 -- Spring, 2016
 *  pgm1:  Linear and nonlinear advection
 *  >>>>> PUT YOUR NAME HERE!  Yes, you!! <<<<<
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include "std2file.h"
using namespace std;

#include <sys/time.h>

#include "prototypes.cuh"
#include "real.h"
#include "dimCube.h"
#include "fatal.hpp"



/*
 * Definitions
 * NX is the number of physical grid points - Not including 'ghost points'
 * bc_width is the number of ghost points - on each side of the physical grid
 * iW is the first C index to use - the left-most physical point
 * iE is the last C index to use - the right-most physical point
 * nxdim is the total array size, including grid points; use for declarations
 * MAXSTEP is the maximum number of time steps to take
 * "c" is [constant] flow speed (e.g. meters/sec) - used only in linear case
 * "dx" is grid spacing: distance (e.g. meters) between grid points in X, Y
 * "name" is a variable containing your name - used to label the plots.
 */



int main(int argc, char *argv[])
{

    double elapsed_time;
    struct timeval tp;
    // Input and output streams
    ifstream infile;
    ofstream outfile;
    // end of Input and output streams

    real dx,dy,dz,dt;

    /* constants for Computer Problem 5 and 6*/
    const real g = 9.81;
    const real thetaBar = 300.;
    
	real deltaTheta[2];      /* two temperature perturbations */
	real deltaV[2];          /* two V perturbations */
	real x0[2],y0[2],z0[2];  /* two locations perturbations */
    real rx[2],ry[2],rz[2];  /* two locations perturbations radius */
	real tstep;

// Arrays and other variables //

    real *ro_u, *ro_w;
    real ***t1, ***t2;
    float *tplot;
    real ***p1, ***p2, ***p3;
    real ***u1, ***u2, ***u3;
    real ***v1, ***v2, ***v3;
    real ***w1, ***w2, ***w3;
	real cs,deltau;
	real ku,kv,kw,kt;
    int nx,ny,nz,iW, iE, jS,jN, kB, kT, nxdim,nydim, nzdim, bc_width=0;
    int nstep,nplot;
    char advection_type;


    /* Parameters and input .................................... */

    if (argc > 1 ) {
        infile.open(argv[1], ios::in);
        if (!infile) {
            cerr <<" Can't open input file " << argv[1] << '!' << endl;
            exit(-1);
        } // end if // 
        std2file(cin,infile);
    } else {
        cout << "Use: " << argv[0] << " filename \n";
        cout << "Example: " << argv[0] << "  data.txt \n"; 
        cout << "advection_scheme can be: L for Lax-Wendroff, C for Crowley, T for Takacs, P for Piecewise, U for Upstream\n";
        exit(-1);
    } // endif //
    
    if (sizeof(real) == 8) {
        cout << "Double precision version\n";
    } else {
        cout << "Single precision version\n";
    } // end if


    cin >> nstep >> nplot >> advection_type;
    cin >> nx >> ny >> nz;
    cin >> dx >> dy >> dz >> dt;
    cin >> cs;
    cin >> ku >> kv >> kw >> kt;
    cin >> deltau;
    cin >> x0[0] >> y0[0] >> z0[0];
    cin >> rx[0] >> ry[0] >> rz[0];
    cin >> deltaTheta[0] >>  deltaV[0];
    cin >> x0[1] >> y0[1] >> z0[1];
    cin >> rx[1] >> ry[1] >> rz[1];
    cin >> deltaTheta[1] >>  deltaV[1];
    


    cout << "Program 6, ATMS 502/CSE 566, Fall 2016\n";
    cout << "Domain (" << setw(3) << nx << ","
                       << setw(3) << ny << ","     
                       << setw(3) << nz << ")"  << endl;

    cout << "dx=" << fixed << setprecision(5) << dx << endl;
    cout << "dy=" << fixed << setprecision(5) << dy << endl;
    cout << "dz=" << fixed << setprecision(5) << dz << endl; 
    cout << "dt=" << fixed << setprecision(5) << dt << endl;
    cout << "cs=" << fixed << setprecision(5) << cs << endl;  
    cout << "KU=" << fixed << setprecision(5) << ku << endl;
    cout << "KV=" << fixed << setprecision(5) << kv << endl;
    cout << "KW=" << fixed << setprecision(5) << kw << endl;
    cout << "KT=" << fixed << setprecision(5) << kt << endl;
    cout << "deltau=" << fixed << setprecision(5) << deltau << endl;
                       
	
    cout << "x_0=" << fixed << setprecision(5) << x0[0] << endl;
    cout << "y_0=" << fixed << setprecision(5) << y0[0] << endl;
    cout << "z_0=" << fixed << setprecision(5) << z0[0] << endl;
    cout << "rx_0=" << fixed << setprecision(5) << rx[0] << endl;
    cout << "ry_0=" << fixed << setprecision(5) << ry[0] << endl;
    cout << "rz_0=" << fixed << setprecision(5) << rz[0] << endl;
    cout << "deltaTheta_0=" << fixed << setprecision(5) << deltaTheta[0] << endl;
    cout << "deltaV_0=" << fixed << setprecision(5) << deltaV[0] << endl;

    cout << "x_1=" << fixed << setprecision(5) << x0[1] << endl;
    cout << "y_1=" << fixed << setprecision(5) << y0[1] << endl;
    cout << "z_1=" << fixed << setprecision(5) << z0[1] << endl;
    cout << "rx_1=" << fixed << setprecision(5) << rx[1] << endl;
    cout << "ry_1=" << fixed << setprecision(5) << ry[1] << endl;
    cout << "rz_1=" << fixed << setprecision(5) << rz[1] << endl;
    cout << "deltaTheta_1=" << fixed << setprecision(5) << deltaTheta[1] << endl;
    cout << "deltaV_1=" << fixed << setprecision(5) << deltaV[1] << endl;

    switch ( advection_type ) {
    case 'l':
    case 'u':
        bc_width = 1;
        break;
    case 'p':
    case 't':
        bc_width = 2;
        break;
    case 'c':
        bc_width = 3;
        break;
    case 'n':
        break;
    default:
        cout << "Integrate: Error, unrecognized advection type " << advection_type << endl;
        exit(-1);
        break;
    } // end switch //
    iW          = bc_width;
    iE          = iW+nx-1;
    jS          = bc_width;
    jN          = jS+ny-1;
    kB          = bc_width;
    kT          = kB+nz-1;
    nxdim       = nx+2*bc_width;
    nydim       = ny+2*bc_width;
    nzdim       = nz+2*bc_width;


    
    
    t1 = dimCube(nzdim,nydim,nxdim);
    t2 = dimCube(nzdim,nydim,nxdim);
	tplot = (float *) malloc( nx * ny * nz * sizeof(float*));

    p1 = dimCube(nz,ny,nx);
    p2 = dimCube(nz,ny,nx);
    p3 = dimCube(nz,ny,nx);


    u1 = dimCube(nzdim,nydim,nxdim+1);
    u2 = dimCube(nzdim,nydim,nxdim+1);
    u3 = dimCube(nzdim,nydim,nxdim+1);
    
    v1 = dimCube(nzdim,nydim+1,nxdim);
    v2 = dimCube(nzdim,nydim+1,nxdim);
    v3 = dimCube(nzdim,nydim+1,nxdim);
    

    w1 = dimCube(nzdim+1,nydim,nxdim);
    w2 = dimCube(nzdim+1,nydim,nxdim);
    w3 = dimCube(nzdim+1,nydim,nxdim);
    

    ro_u = (real *) malloc(nzdim*sizeof(real));
    ro_w = (real *) malloc((nzdim+1)*sizeof(real));
    
    
    // end of allocating memory for arrays and matrices


    gettimeofday(&tp,NULL);
    elapsed_time = -(tp.tv_sec*1.0e6 + tp.tv_usec);  

    /*
     * Set and plot the initial condition
     */
     
    ic(t1,p1,u1,v1,w1,ro_u,ro_w,dx,dy,dz,deltau,iW,iE,jS,jN,kB,kT,bc_width,nx,ny,nz,x0,y0,z0,deltaTheta,deltaV,thetaBar,rx,ry,rz,g);
    
    bc(u1,v1,w1,iW,iE,jS,jN,kB,kT,bc_width);
    bc4T(t1,t2,iW,iE,jS,jN,kB,kT,bc_width);
    // acctually only the BC of t1 should be copyed to t2 
    
    for (int level=0; level<nzdim; ++level) {
        for (int row=0; row<nydim; ++row) {
            for (int col=0; col<nxdim; ++col) {
                t2[level][row][col]=t1[level][row][col];
            } // end for //
        } // end for //
    } // end for //
    
    // Ploting initial conditions //
    
    for (int col=iW; col<=iE; ++col) {
        for (int row=jS; row<=jN; ++row) {
            for (int level=kB; level<=kT; ++level) {
                tplot[(col-iW)*nz*ny + (row-jS)*nz + (level-kB)] = (float) (t1[level][row][col] - thetaBar);
            } // end for //
        } // end for //
    } // end for // 
   
    putfield((char *) "T",0.0,tplot, nx,ny,nz);
    
    
    for (int col=iW; col<=iE; ++col) {
        for (int row=jS; row<=jN; ++row) {
            for (int level=kB; level<=kT; ++level) {
                tplot[(col-iW)*nz*ny + (row-jS)*nz + (level-kB)] = (float) (0.5*(u1[level][row][col]+u1[level][row][col+1]));
            } // end for //
        } // end for //
    } // end for // 
   
    putfield((char *) "U",0.0,tplot, nx,ny,nz);
    
    
    for (int col=iW; col<=iE; ++col) {
        for (int row=jS; row<=jN; ++row) {
            for (int level=kB; level<=kT; ++level) {
                tplot[(col-iW)*nz*ny + (row-jS)*nz + (level-kB)] = (float) (0.5*(v1[level][row][col]+v1[level][row+1][col]));
            } // end for //
        } // end for //
    } // end for //     
    
    putfield((char *) "V",0.0,tplot, nx,ny,nz);
    
    
    for (int col=iW; col<=iE; ++col) {
        for (int row=jS; row<=jN; ++row) {
            for (int level=kB; level<=kT; ++level) {
                tplot[(col-iW)*nz*ny + (row-jS)*nz + (level-kB)] = (float) (0.5*(w1[level][row][col]+w1[level+1][row][col]));
            } // end for //
        } // end for //
    } // end for //
     
    putfield((char *) "W",0.0,tplot, nx,ny,nz);

   
    for (int col=0; col<nx; ++col) {
        for (int row=0; row<ny; ++row) {
            for (int level=0; level<nz; ++level) {
                tplot[col*nz*ny + row*nz + level] = (float) (p1[level][row][col]);
            } // end for //
        } // end for //
    } // end for //
    
    putfield((char *) "P",0.0,tplot, nx,ny,nz);

    // end of Plotting initial conditions //


    /*
     * .. Integrate .....................................................
     */
    
     
    tstep=dt;     
    for (int n=1 ; n<=nstep; n++) {
        
        for (int level=0; level<nzdim; ++level) {
            for (int row=0; row<nydim; ++row) {
                for (int col=0; col<=nxdim; ++col) {
                    u3[level][row][col]=u1[level][row][col];
                } // end for //
            } // end for //
        } // end for //

        
        for (int level=0; level<nzdim; ++level) {
            for (int row=0; row<=nydim; ++row) {
                for (int col=0; col<nxdim; ++col) {
                    v3[level][row][col]=v1[level][row][col];
                } // end for //
            } // end for //
        } // end for //
        
        
        for (int level=0; level<=nzdim; ++level) {
            for (int row=0; row<nydim; ++row) {
                for (int col=0; col<nxdim; ++col) {
                    w3[level][row][col]=w1[level][row][col];
                } // end for //
            } // end for //
        } // end for //
        
        //  . . . Compute values at next step                           //


        advection(t2,u3,v3,w3,t1,u2,v2,w2,dx,dy,dz,dt,iW,iE,jS,jN,kB,kT,nxdim,nydim,nzdim,tstep,advection_type);

        diffusion(t2,u3,v3,w3,t1,u1,v1,w1,dx*dx,dy*dy,dz*dz,iW,iE,jS,jN,kB,kT,dt,tstep,ku,kv,kw,kt);// KU,KV,KW,KT);

        pgf1(u3,v3,w3,t2,p1,ro_u,ro_w,dx,dy,dz,iW,iE,jS,jN,kB,kT,tstep, g,thetaBar,bc_width);
        bc(u3,v3,w3,iW,iE,jS,jN,kB,kT,bc_width);
        pgf2(p3,u3,v3,w3,t2,p1,ro_u,ro_w,dx,dy,dz,iW,iE,jS,jN,kB,kT,tstep,cs*cs,bc_width);

        //  . . . Set boundary conditions for T                     //
        bc4T(t2,t1,iW,iE,jS,jN,kB,kT,bc_width);


        if (n%nplot == 0) {
        
            
            for (int col=iW; col<=iE; ++col) {
                for (int row=jS; row<=jN; ++row) {
                    for (int level=kB; level<=kT; ++level) {
                        tplot[(col-iW)*nz*ny + (row-jS)*nz + (level-kB)] = (float) (t2[level][row][col] - thetaBar);
                    } // end for //
                } // end for //
            } // end for // 
            
            putfield((char *) "T",(float) dt*n,tplot, nx,ny,nz);
            
            
            for (int col=iW; col<=iE; ++col) {
                for (int row=jS; row<=jN; ++row) {
                    for (int level=kB; level<=kT; ++level) {
                        tplot[(col-iW)*nz*ny + (row-jS)*nz + (level-kB)] = (float) (0.5*(u3[level][row][col]+u3[level][row][col+1]));
                    } // end for //
                } // end for //
            } // end for // 
             
            putfield((char *) "U",(float) dt*n,tplot, nx,ny,nz);

            
            for (int col=iW; col<=iE; ++col) {
                for (int row=jS; row<=jN; ++row) {
                    for (int level=kB; level<=kT; ++level) {
                        tplot[(col-iW)*nz*ny + (row-jS)*nz + (level-kB)] = (float) (0.5*(v3[level][row][col]+v3[level][row+1][col]));
                    } // end for //
                } // end for //
            } // end for //     
            
            putfield((char *) "V",(float) dt*n,tplot, nx,ny,nz);
            
            
            for (int col=iW; col<=iE; ++col) {
                for (int row=jS; row<=jN; ++row) {
                    for (int level=kB; level<=kT; ++level) {
                        tplot[(col-iW)*nz*ny + (row-jS)*nz + (level-kB)] = (float) (0.5*(w3[level][row][col]+w3[level+1][row][col]));
                    } // end for //
                } // end for //
            } // end for //
            
            putfield((char *) "W",(float) dt*n,tplot, nx,ny,nz);

            
            for (int col=0; col<nx; ++col) {
                for (int row=0; row<ny; ++row) {
                    for (int level=0; level<nz; ++level) {
                        tplot[col*nz*ny + row*nz + level] = (float) (p3[level][row][col]);
                    } // end for //
                } // end for //
            } // end for //
            
            putfield((char *) "P",(float) dt*n,tplot, nx,ny,nz);
        
        } // end if //

        //  . . . Do array update at end of time step                   //
        if (n == 1) { 
            tstep=2.0*dt;
            update(&p2,&p3);
            update(&u2,&u3);
            update(&v2,&v3);
            update(&w2,&w3);
            
        } else {
            update(&p1,&p2);
            update(&p2,&p3);
            
            update(&u1,&u2);
            update(&u2,&u3);
            
            update(&v1,&v2);
            update(&v2,&v3);
            
            update(&w1,&w2);
            update(&w2,&w3);
        } // end if //
        update(&t1,&t2);

    }	// end of time loop n = 1,...,nstep //
        
    gettimeofday(&tp,NULL);
    elapsed_time += (tp.tv_sec*1.0e6 + tp.tv_usec);
    cout << "\n\nIt tooks " << setw(14) << scientific << setprecision(6) 
         << elapsed_time*1.0e-6 << " seconds to finish\n";

    if (sizeof(real) == 8) {
        cout << "Double precision version\n";
    } else {
        cout << "Single precision version\n";
    } // end if


    
    freeCube(&t2);
    freeCube(&t1);

    freeCube(&p1);
    freeCube(&p2);
    freeCube(&p3);
    
    freeCube(&u1);
    freeCube(&u2);
    freeCube(&u3);

	freeCube(&v1);
	freeCube(&v2);
	freeCube(&v3);

    freeCube(&w1);
    freeCube(&w2);
    freeCube(&w3);

	free(tplot);
	free(ro_u);
	free(ro_w);

    
    return 0;
} // end main() //

