#ifndef REAL
    #define FLOAT
    //#define DOUBLE

    #ifdef FLOAT
       typedef float real;
    #else
       typedef double real;
    #endif
    #define REAL
#endif

