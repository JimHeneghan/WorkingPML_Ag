#ifndef _P_MACRO_H
#define _P_MACRO_H
#include <complex.h>

#define SizeX0 10
#define SizeY0 10
/****************************************************/
/*********************IMPORTANT!!********************/
/*SizeZ0 must include PML layers AND Simulationspace*/
#define SizeZ0 830
#define PMLs 115
/****************************************************/
#define MT 4e4

#define DDx0 10e-9
#define DDy0 10e-9
#define DDz0 1e-9

#define pPMLHx_1(PP)       pPMLhx_1[PP]
#define pPMLHx_2(PP)       pPMLhx_2[PP]

#define pPMLHz_2(PP)       pPMLhz_2[PP]
#define pPMLHz_3(PP)       pPMLhz_3[PP]

#define pPMLDx_1(PP)       pPMLdx_1[PP]
#define pPMLDx_2(PP)       pPMLdx_2[PP]

#define pPMLDz_2(PP)       pPMLdz_2[PP]
#define pPMLDz_3(PP)       pPMLdz_3[PP]

#define pRXSensor(MM)      prxsensor[MM]
#define pIXSensor(MM)      pixsensor[MM]
#define pTXSensor(MM)      ptxsensor[MM]

#define pRYSensor(MM)      prysensor[MM]
#define pIYSensor(MM)      piysensor[MM]
#define pTYSensor(MM)      ptysensor[MM]

#define pRZSensor(MM)      przsensor[MM]
#define pIZSensor(MM)      pizsensor[MM]
#define pTZSensor(MM)      ptzsensor[MM]

#define pHx(MM,NN,PP)      phx[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pHy(MM,NN,PP)      phy[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pHz(MM,NN,PP)      phz[((MM)*SizeY0 + NN)*SizeZ0 + PP]

#define pDx(MM,NN,PP)      pdx[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pDy(MM,NN,PP)      pdy[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pDz(MM,NN,PP)      pdz[((MM)*SizeY0 + NN)*SizeZ0 + PP]

#define pEx(MM,NN,PP)      pex[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pEy(MM,NN,PP)      pey[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pEz(MM,NN,PP)      pez[((MM)*SizeY0 + NN)*SizeZ0 + PP]

#define pCurlHx(MM,NN,PP)      pCurlhx[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pCurlHy(MM,NN,PP)      pCurlhy[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pCurlHz(MM,NN,PP)      pCurlhz[((MM)*SizeY0 + NN)*SizeZ0 + PP]

#define pCurlEx(MM,NN,PP)  pCurlex[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pCurlEy(MM,NN,PP)  pCurley[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pCurlEz(MM,NN,PP)  pCurlez[((MM)*SizeY0 + NN)*SizeZ0 + PP]

#define pICEz(MM,NN,PP)      pICez[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pICHz(MM,NN,PP)      pIChz[((MM)*SizeY0 + NN)*SizeZ0 + PP]

#define pSxx(MM,NN,PP)      psxx[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pSyy(MM,NN,PP)      psyy[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pSz(MM,NN,PP)       psz[((MM)*SizeY0 + NN)*SizeZ0 + PP]

#define pS1xx(MM,NN,PP)     ps1xx[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pS1yy(MM,NN,PP)     ps1yy[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pS1z(MM,NN,PP)      ps1z[((MM)*SizeY0 + NN)*SizeZ0 + PP]

#define pS2xx(MM,NN,PP)     ps2xx[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pS2yy(MM,NN,PP)     ps2yy[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pS2z(MM,NN,PP)      ps2z[((MM)*SizeY0 + NN)*SizeZ0 + PP]

#define pE1x(MM,NN,PP)      pe1x[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pE1y(MM,NN,PP)      pe1y[((MM)*SizeY0 + NN)*SizeZ0 + PP]
#define pE1z(MM,NN,PP)      pe1z[((MM)*SizeY0 + NN)*SizeZ0 + PP]

#define pMedia(MM,NN,PP)    pmedia[((MM)*SizeY0 + NN)*(SizeZ0) + PP]

#endif
