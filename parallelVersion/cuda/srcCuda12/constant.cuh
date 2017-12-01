#ifdef MAIN

    __constant__ int nx,ny,nz,iW, iE, jS,jN, kB, kT, nxdim,nydim, nzdim, bcw;
    
    __constant__ char advection_type;

    __constant__ real KU,KV,KW,KTemp,g, thetaBar,dx,dy,dz,dt,cs2,tstep;
    
    __constant__ real *ro_u, *ro_w;

#else

    extern __constant__ int nx,ny,nz,iW, iE, jS,jN, kB, kT, nxdim,nydim, nzdim, bcw;

    extern __constant__ char advection_type;

    extern __constant__ real KU,KV,KW,KTemp, g, thetaBar,dx,dy,dz,dt,cs2,tstep;

    extern __constant__ real *ro_u, *ro_w;

#endif

#define TILESIZEX 16
#define TILESIZEY 16

