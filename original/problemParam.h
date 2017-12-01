#include <float.h>
#define TEST_A3
/*

#define OFFICIAL
#define TEST_A1
#define TEST_A2
#define TEST_A3
#define TEST_B
#define TEST_C
#define TEST_D1
#define TEST_D2
#define TEST_E
#define TEST_F
*/


#if defined OFFICIAL //_HIGH
    #define NX 300
    #define NY 300   
    #define NZ 75
    #define DX 50
    #define DY 50
    #define DZ 50
    #define CS 60.0
    #define KU 40.0
    #define KV 40.0
    #define KW 40.0
    #define KT 5.0        
    #define DT 0.15
    #define DELTAU 2.0
    #define DELTAV_0 -40.0
    #define DELTAV_1  40.0    
    /* 1st perturbation case */
    #define DELTATHETA_0 -25.0
    #define X0_0 25
    #define Y0_0 7525
    #define Z0_0 1525
    #define XRADIUS_0 3500
    #define YRADIUS_0 FLT_MAX
    #define ZRADIUS_0 1750
    /* 2nd perturbation case */
    #define DELTATHETA_1 -25.0
    #define X0_1 14975
    #define Y0_1 7525
    #define Z0_1 1525
    #define XRADIUS_1 3500
    #define YRADIUS_1 FLT_MAX
    #define ZRADIUS_1 1750
#elif defined OFFICIAL_LOW
    #define NX 151
    #define NY 151   
    #define NZ 35
    #define DX 150
    #define DY 150  
    #define DZ 150
    #define CS 60.0
    #define KU 50.0
    #define KV 50.0
    #define KW 50.0
    #define KT 5.0        
    #define DT 0.50
    #define DELTAU 2.0
    #define DELTAV_0 -35.0
    #define DELTAV_1  35.0    
    /* 1st perturbation case */
    #define DELTATHETA_0 -25.0
    #define X0_0 75
    #define Y0_0 11325
    #define Z0_0 1575
    #define XRADIUS_0 4000
    #define YRADIUS_0 FLT_MAX
    #define ZRADIUS_0 2000
    /* 2nd perturbation case */
    #define DELTATHETA_1 -25.0
    #define X0_1 22575
    #define Y0_1 11325
    #define Z0_1 1575
    #define XRADIUS_1 4000
    #define YRADIUS_1 FLT_MAX
    #define ZRADIUS_1 2000
#elif defined OFFICIAL_LOW2
    #define NX 151
    #define NY 151   
    #define NZ 35
    #define DX 37.5
    #define DY 37.5  
    #define DZ 37.5
    #define CS 60.0
    #define KU 40.0
    #define KV 40.0
    #define KW 40.0
    #define KT 5.0        
    #define DT 0.10
    #define DELTAU 2.0
    /* 1st perturbation case */
    #define X0_0 18.75
    #define Y0_0 2831.25
    #define Z0_0 393.75 
    #define DELTATHETA_0 -20.0
    #define DELTAV_0 -30.0
    #define XRADIUS_0 750
    #define YRADIUS_0 FLT_MAX
    #define ZRADIUS_0 750
    /* 2nd perturbation case */
    #define X0_1 5643.75
    #define Y0_1 2831.25
    #define Z0_1 393.75 
    #define DELTATHETA_1 -20.0
    #define DELTAV_1 30.0
    #define XRADIUS_1 750
    #define YRADIUS_1 FLT_MAX
    #define ZRADIUS_1 750
#elif defined TEST_A1
    #define NX 33
    #define NY 33    
    #define NZ 16
    #define DX 500
    #define DY 500   
    #define DZ 500
    #define DT 1.00
    #define CS 100.0
    /* 1st perturbation case */
    #define X0_0 8250
    #define Y0_0 8250
    #define Z0_0 4250
    #define DELTATHETA_0 -20.0
    #define XRADIUS_0 4000
    #define YRADIUS_0 4000
    #define ZRADIUS_0 4000
    /* 2nd perturbation case */
    #define X0_1 8250
    #define Y0_1 8250
    #define Z0_1 4250
    #define DELTATHETA_1 0.0
    #define XRADIUS_1 4000
    #define YRADIUS_1 4000
    #define ZRADIUS_1 4000
    
    #define KU 00.0
    #define KV 00.0
    #define KW 00.0
    #define KT 50.0    
    
#elif defined TEST_A2
    #define NX 33
    #define NY 33    
    #define NZ 16
    #define DX 500
    #define DY 500   
    #define DZ 500
    #define DT 1.00
    #define CS 100.0
    /* 1st perturbation case */
    #define X0_0 8250
    #define Y0_0 8250
    #define Z0_0 4250
    #define DELTATHETA_0 -20.0
    #define XRADIUS_0 4000
    #define YRADIUS_0 4000
    #define ZRADIUS_0 4000
    /* 2nd perturbation case */
    #define X0_1 8250
    #define Y0_1 8250
    #define Z0_1 4250
    #define DELTATHETA_1 0.0
    #define XRADIUS_1 4000
    #define YRADIUS_1 4000
    #define ZRADIUS_1 4000
    
    #define KU 00.0
    #define KV 00.0
    #define KW 00.0
    #define KT 50.0    
    
#elif defined TEST_A3
    #define NX 33
    #define NY 33    
    #define NZ 16
    #define DX 500
    #define DY 500   
    #define DZ 500
    #define DT 1.00
    #define CS 100.0
    /* 1st perturbation case */
    #define X0_0 8250
    #define Y0_0 8250
    #define Z0_0 4250
    #define DELTATHETA_0 -20.0
    #define XRADIUS_0 4000
    #define YRADIUS_0 4000
    #define ZRADIUS_0 4000
    /* 2nd perturbation case */
    #define X0_1 8250
    #define Y0_1 8250
    #define Z0_1 4250
    #define DELTATHETA_1 0.0    
    #define XRADIUS_1 4000
    #define YRADIUS_1 4000
    #define ZRADIUS_1 4000

    #define KU 00.0
    #define KV 00.0
    #define KW 00.0
    #define KT 50.0    

#elif defined TEST_B
    #define NX 33
    #define NY 33    
    #define NZ 16
    #define DX 500
    #define DY 500   
    #define DZ 500
    #define DT 1.00
    #define CS 100.0
    /* 1st perturbation case */
    #define X0_0 8250
    #define Y0_0 8250
    #define Z0_0 4250
    #define DELTATHETA_0 -20.0
    #define XRADIUS_0 4000
    #define YRADIUS_0 4000
    #define ZRADIUS_0 4000
    /* 2nd perturbation case */
    #define X0_1 8250
    #define Y0_1 8250
    #define Z0_1 4250
    #define DELTATHETA_1 0.0    
    #define XRADIUS_1 4000
    #define YRADIUS_1 4000
    #define ZRADIUS_1 4000

    #define KU 00.0
    #define KV 00.0
    #define KW 00.0
    #define KT 50.0    

#elif defined TEST_C
    #define NX 33
    #define NY 33    
    #define NZ 16
    #define DX 500
    #define DY 500   
    #define DZ 500
    #define DT 1.00
    #define CS 100.0
    /* 1st perturbation case */
    #define X0_0 8250
    #define Y0_0 8250
    #define Z0_0 4250
    #define DELTATHETA_0 -20.0
    #define XRADIUS_0 4000
    #define YRADIUS_0 4000
    #define ZRADIUS_0 4000
    /* 2nd perturbation case */
    #define X0_1 8250
    #define Y0_1 8250
    #define Z0_1 4250
    #define DELTATHETA_1 0.0    
    #define XRADIUS_1 4000
    #define YRADIUS_1 4000
    #define ZRADIUS_1 4000
    
    #define KU 00.0
    #define KV 00.0
    #define KW 00.0
    #define KT 50.0
    
#elif (defined TEST_D1 || defined TEST_D2) 
    #define NX 33
    #define NY 33    
    #define NZ 16
    #define DX 500
    #define DY 500   
    #define DZ 500
    #define DT 1.00
    #define CS 100.0
    /* 1st perturbation case */
    #define X0_0 8250
    #define Y0_0 8250
    #define Z0_0 4250
    #define DELTATHETA_0 -20.0
    #define XRADIUS_0 4000
    #define YRADIUS_0 4000
    #define ZRADIUS_0 4000
    /* 2nd perturbation case */
    #define X0_1 8250
    #define Y0_1 8250
    #define Z0_1 4250
    #define DELTATHETA_1 0.0   
    #define XRADIUS_1 4000
    #define YRADIUS_1 4000
    #define ZRADIUS_1 4000
    
    #define KU 50.0
    #define KV 50.0
    #define KW 50.0
    #define KT  5.0
#elif (defined TEST_E) 
    #define NX 33
    #define NY 33    
    #define NZ 16
    #define DX 500
    #define DY 500   
    #define DZ 500
    #define DT 1.00
    #define CS 100.0
    /* 1st perturbation case */
    #define X0_0 250
    #define Y0_0 8250
    #define Z0_0 4250
    #define DELTATHETA_0 -20.0
    #define DELTAV_0 -25.0
    #define XRADIUS_0 4000
    #define YRADIUS_0 FLT_MAX
    #define ZRADIUS_0 4000
    /* 2nd perturbation case */
    #define X0_1 16250
    #define Y0_1 8250
    #define Z0_1 4250
    #define DELTATHETA_1 -20.0
    #define DELTAV_1 25.0    
    #define XRADIUS_1 4000
    #define YRADIUS_1 FLT_MAX
    #define ZRADIUS_1 4000
    
    #define KU 50.0
    #define KV 50.0
    #define KW 50.0
    #define KT  5.0
#elif (defined TEST_F) 
    #define NX 53
    #define NY 53    
    #define NZ 16
    #define DX 500
    #define DY 500   
    #define DZ 500
    #define DT 1.00
    #define CS 100.0
    #define DELTAU 10.0    
    /* 1st perturbation case */
    #define X0_0 250
    #define Y0_0 10750
    #define Z0_0 1250
    #define DELTATHETA_0 -10.0
    #define DELTAV_0 -25.0
    #define XRADIUS_0 4000
    #define YRADIUS_0 FLT_MAX
    #define ZRADIUS_0 4000
    /* 2nd perturbation case */
    #define X0_1 16250
    #define Y0_1 8250
    #define Z0_1 4250
    #define DELTATHETA_1 0.0
    #define DELTAV_1 0.0    
    #define XRADIUS_1 4000
    #define YRADIUS_1 FLT_MAX
    #define ZRADIUS_1 4000
    
    #define KU 35.0
    #define KV 35.0
    #define KW 35.0
    #define KT  2.5
#endif

