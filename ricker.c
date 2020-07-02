#include "ezinc.h"
#include "p-macro.h"
#include <openacc.h>

static double cdtds, ppw = 0;

/* initialize source-function variables */
void ezIncInit(Grid *g){

  ppw = 30;
  cdtds = Cdtds;
  return;
}

/* calculate source function at given time and location */
void TFSF_E(double time, Grid *g) {
  double arg;
  double rick, sig, dt, ddx, t0, tau, Gauss, c0, ExSr, Eysr, Nlam;
  int mm, nn, pp;

  int tc = SizeX0*SizeY0*SizeZ0;
  complex double *pCurlex = g->Curlex;
  complex double *pCurley = g->Curley;
  complex double *pCurlez = g->Curlez;
  complex double *pex = g->ex;
  complex double *pey = g->ey;
  complex double *pez = g->ez;

  double DDx = DDx0;
  double DDy = DDy0;
  double DDz = DDz0;

  if (ppw <= 0) {
    fprintf(stderr,
       "ezInc: ezIncInit() must be called before ezInc.\n"
       "       Points per wavelength must be positive.\n");
    exit(-1);
  }
  c0 = 3e8;
  ddx = DDz;
  // dt = ddx/(2*c0);
  // //sig = 10;
  // tau = 1000*dt;
  // t0 = 2*tau;
  Nlam = 250.0;
  // omega = 90e12*2.0*M_PI;
  // arg = ((time*dt - t0)/(tau)) ;//M_PI * ((cdtds * time - location) / ppw - 1.0);
  // arg = arg * arg;

  // C1 =  cos(omega*(time*dt - t0));
  //printf("ricker is %g \n", rick);
  //return (exp(-arg));//(1.0 - 2.0 * arg) * exp(-arg);

  // Gauss = C1*(exp(-arg));
  arg = (time*M_PI)/(3*Nlam) - 10.0;
  arg = arg*arg;
  Gauss = exp(-0.5*arg)*cos((time*M_PI/Nlam) - 30.0);

  pp = 15 + PMLs;

  // #pragma acc enter data copyin(pex[0:tc], pey[0:tc], pez[0:tc], pCurlex[0:tc], pCurley[0:tc], pCurlez[0:tc])
  // {
    #pragma acc parallel loop async(1) independent collapse(2) present(pCurlex[0:tc], pey[0:tc],pez[0:tc])
    for (mm = 0; mm < SizeX0; mm++){
      for (nn = 0; nn < SizeY0; nn++){
        pCurlEx(mm, nn, pp) = pCurlEx(mm, nn, pp) + Gauss/DDz;
      }
    }

    // #pragma acc parallel loop async(2) independent present(pCurlex[0:tc], pey[0:tc],pez[0:tc])
    // for (mm = 0; mm < SizeX0; mm++){
    //  pCurlEx(mm, SizeY0 - 1, pp) = (pEz(mm, 0, pp) - pEz(mm, SizeY0 - 1, pp))/DDy
    //   - (pEy(mm, SizeY0 - 1, pp + 1) - pEy(mm, SizeY0 -1, pp))/DDz + Gauss/DDz;
    // }

    #pragma acc parallel loop async(3) independent collapse(2) present(pCurley[0:tc], pex[0:tc],pez[0:tc])
    for (mm = 0; mm < SizeX0; mm++){
      for (nn = 0; nn < SizeY0; nn++){
        pCurlEy(mm, nn, pp) = pCurlEy(mm, nn, pp) - Gauss/DDz;
      }
    }  
    // #pragma acc parallel loop async(4) independent present(pCurley[0:tc], pex[0:tc],pez[0:tc])
    // for (nn = 0; nn < SizeY0; nn++){
    //     pCurlEy(SizeX0 - 1, nn, pp) = (pEx(SizeX0 - 1, nn, pp + 1) - pEx(SizeX0 - 1, nn, pp))/DDz
    //      - (pEz(0, nn, pp) - pEz(SizeX0 - 1, nn, pp))/DDx- Gauss/DDz;
    // }
    // #pragma acc exit data copyout(pex[0:tc], pey[0:tc], pez[0:tc], pCurlex[0:tc], pCurley[0:tc], pCurlez[0:tc])
  // }
    #pragma acc wait
    return;
}


void TFSF_H(double time, Grid *g) {
  double arg;
  double rick, sig, dt, ddx, t0, tau, Gauss, c0, del_t, Hysr, Hxsr, Nlam;
  int mm, nn, pp;

  double imp0 = 377.0;
  int tc = SizeX0*SizeY0*SizeZ0;
  complex double *phx = g->hx;
  complex double *phy = g->hy;
  complex double *phz = g->hz;
  
  complex double *pCurlhx = g->Curlhx;
  complex double *pCurlhy = g->Curlhy;

  double DDx = DDx0;
  double DDy = DDy0;
  double DDz = DDz0;

  if (ppw <= 0) {
    fprintf(stderr,
       "ezInc: ezIncInit() must be called before ezInc.\n"
       "       Points per wavelength must be positive.\n");
    exit(-1);
  }
  c0 = 3e8;
  ddx = DDz;
  // dt = ddx/(2*c0);
  // tau= 1000*dt;
  // t0 = 2*tau;
  // omega = 90e12*2.0*M_PI;
  /*assuming cubic Yee cells*/
  //del_t = -1*(ddx/2*c0) + (dt/2);
  //del_t = -1 + 0.5;
  del_t = -0.5*dt;
  Nlam = 250.0;

  // arg = ((time - t0 + del_t)/(tau)) ;//M_PI * ((cdtds * time - location) / ppw - 1.0);
  // arg = arg * arg;

  // C1 =  cos(omega*(time*dt - t0 + del_t));

  // Gauss = C1*(exp(-arg));

  arg = ((time - 0.5)*M_PI)/(3*Nlam) - 10.0;
  arg = arg*arg;
  Gauss = exp(-0.5*arg)*cos(((time - 0.5)*M_PI/Nlam) - 30.0);

  pp = 15 + PMLs;
  // Hxsr = -1*Gauss;
  // Hysr =    Gauss;
  //printf("0 CurlHx(20, 20, 30) is %g \n", CurlHx(20, 20, 30));

  // #pragma acc enter data copyin(pCurlhx[0:tc], pCurlhy[0:tc])
  // {
  #pragma acc parallel loop async(1) independent collapse(2) present(pCurlhx[0:tc],phz[0:tc], phy[0:tc])  
    for (mm = 0; mm < SizeX0; mm++){
      for (nn = 0; nn < SizeY0; nn++){
      pCurlHx(mm, nn, pp) = pCurlHx(mm, nn, pp) + Gauss/DDz;
      }
    }
  // #pragma acc parallel loop async(2) independent present(pCurlhx[0:tc],phz[0:tc], phy[0:tc])
  //   for (mm = 0; mm < SizeX0; mm++){
  //       pCurlHx(mm, 0, pp) = pCurlHx(mm, nn, pp) + Gauss/DDz;
  //   }

    /*for periodic boundary may need to add in the TFSF term explicitly*/
  #pragma acc parallel loop async(3) independent collapse(2) present(pCurlhy[0:tc],phx[0:tc], phz[0:tc])
    for (mm = 0; mm < SizeX0; mm++){
        for (nn = 0; nn < SizeY0; nn++){
          pCurlHy(mm, nn, pp) = pCurlHy(mm, nn, pp)  + Gauss/DDz;
        }
    }
  // #pragma acc parallel loop async(4) independent present(pCurlhy[0:tc],phx[0:tc], phz[0:tc])
  //   for (nn = 0; nn < SizeY0; nn++){
  //       pCurlHy(0, nn, pp) = pCurlHy(mm, nn, pp)  + Gauss/DDz;
  //   }
  // #pragma acc exit data copyout(phx[0:tc],phy[0:tc],phz[0:tc],pCurlhx[0:tc], pCurlhy[0:tc], pCurlhz[0:tc])

    //printf("CurlHx(20, 20, 30) is %g \n \n", CurlHx(20, 20, 30));
  #pragma acc wait
    return;
}
