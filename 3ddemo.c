/* 3D simulation with dipole source at center of grid. */

#include "fdtd-alloc.h"
#include "fdtd-macro.h"
#include "p-macro.h" 
#include "fdtd-proto.h"
#include "ezinc.h"
#include <stdio.h>
#include <complex.h>
#include <openacc.h>



int main() {
  Grid *g;
  
  ALLOC_1D(g, 1, Grid); // allocate memory for grid structure

  double delta1, delta2;

  sensorInit(g);
  
  gridInit(g);        // initialize 3D grid
  
  ezIncInit(g); 
  // snapshot3dInit(g);  // initialize snapshots
  

  printf("Init done \n");

  
  int tc = SizeX0*SizeY0*SizeZ0;

  /* Redefining pointers to arrays to use in the GPU calculations so that they are not part of a struct */
  complex double *phx = g->hx;
  complex double *phy = g->hy;
  complex double *phz = g->hz;

  complex double *pex = g->ex;
  complex double *pey = g->ey;
  complex double *pez = g->ez;

  complex double *pdx = g->dx;
  complex double *pdy = g->dy;
  complex double *pdz = g->dz;

  complex double *pICez = g->ICez;
  complex double *pIChz = g->IChz;

  complex double *pPMLhx_1 = g->PMLhx_1;
  complex double *pPMLhx_2 = g->PMLhx_2;

  complex double *pPMLhz_2 = g->PMLhz_2;
  complex double *pPMLhz_3 = g->PMLhz_3;

  complex double *pPMLdx_1 = g->PMLdx_1;
  complex double *pPMLdx_2 = g->PMLdx_2;

  complex double *pPMLdz_2 = g->PMLdz_2;
  complex double *pPMLdz_3 = g->PMLdz_3;

  complex double *pCurlhx = g->Curlhx;
  complex double *pCurlhy = g->Curlhy;
  complex double *pCurlhz = g->Curlhz;

  complex double *pCurlex = g->Curlex;
  complex double *pCurley = g->Curley;
  complex double *pCurlez = g->Curlez;

  complex double *psxx = g->sxx;
  complex double *psyy = g->syy;
  complex double *psz  = g->sz;

  complex double *ps1xx = g->s1xx;
  complex double *ps1yy = g->s1yy;
  complex double *ps1z  = g->s1z; 
  
  complex double *ps2xx = g->s2xx;
  complex double *ps2yy = g->s2yy;
  complex double  *ps2z = g->s2z;

  complex double *pe1x = g->e1x;
  complex double *pe1y = g->e1y;
  complex double *pe1z = g->e1z;

  complex double *prxsensor = g->rxsensor;
  complex double *pixsensor = g->ixsensor;
  complex double *ptxsensor = g->txsensor;

  complex double *prysensor = g->rysensor;
  complex double *piysensor = g->iysensor;
  complex double *ptysensor = g->tysensor;

  complex double *przsensor = g->rzsensor;
  complex double *pizsensor = g->izsensor;
  complex double *ptzsensor = g->tzsensor;

  int *pmedia = g->media;

/* PML arrays totalling 8x16 bytes of data per Z dimension cell***************/
/* Could be rewritten to only be needed for the actual number of PML layers***/
/* reducing the data required to 8x16 bytes of data per PML layer*************/
#pragma acc enter data copyin(pPMLhx_1[0:SizeZ0],pPMLhx_2[0:SizeZ0],pPMLhz_2[0:SizeZ0],pPMLhz_3[0:SizeZ0])
#pragma acc enter data copyin(pPMLdx_1[0:SizeZ0],pPMLdx_2[0:SizeZ0],pPMLdz_2[0:SizeZ0],pPMLdz_3[0:SizeZ0])

/* Integration terms for the PML update equations*****************************/
/* Currenttly these require 2x16 bytes of data total Yee cell count***********/
/* Could be rewritten to only be needed for the actual number of PML layers***/
/* reducing the data required to 2x16 bytes of data per XxYxPML layer*********/
#pragma acc enter data copyin(pICez[:tc], pIChz[:tc])

/* Field arrays totalling 9x16 bytes per total Yee cell count ****************/
#pragma acc enter data copyin(pex[0:tc],pey[0:tc],pez[0:tc])
#pragma acc enter data copyin(phx[0:tc],phy[0:tc],phz[0:tc])
#pragma acc enter data copyin(pdx[0:tc],pdy[0:tc],pdz[0:tc])

/* Curl arrays totalling 6x16 bytes per total Yee cell count *****************/
/* Possible to delete and create these during the simulation******************/
/* Reducing the memory requierments to only one of the below arrays at a time*/
#pragma acc enter data copyin(pCurlex[0:tc],pCurley[0:tc],pCurlez[0:tc])
#pragma acc enter data copyin(pCurlhx[0:tc],pCurlhy[0:tc],pCurlhz[0:tc])

/* Lorentz Z Transform update arrays for diagonal permittivity****************/
/* arrays totalling 12x16 bytes per total Yee cell count *********************/
/* as written requiers */
#pragma acc enter data copyin(pe1x[:tc], pe1y[:tc], pe1z[:tc])
#pragma acc enter data copyin(psxx[:tc], psyy[:tc], psz[:tc])
#pragma acc enter data copyin(ps1xx[:tc], ps1yy[:tc], ps1z[:tc])
#pragma acc enter data copyin(ps2xx[:tc], ps2yy[:tc], ps2z[:tc])


/* Media characterization array to tell the E field update which set of ******/
/* update equations to use. Currenttly costs 4 bytes per Yee cell count ******/
#pragma acc enter data copyin(pmedia[:tc])


/* Time domain sensor arrays. currenttly cost 9x16x(no of time steps) of data */
/* Could be used every (int) number of steps to reduce memory requierments ***/
/* Could be replaced by a direct FT function that would cost *****************/
/* 16 bytes/sensor/frequency *************************************************/
/* or SizeXxSizeYx16 bytes per sensor if used over an entire plane per freq **/
#pragma acc enter data copyin(prxsensor[0:MT], pixsensor[0:MT], ptxsensor[0:MT])
#pragma acc enter data copyin(prysensor[0:MT], piysensor[0:MT], ptysensor[0:MT])
#pragma acc enter data copyin(przsensor[0:MT], pizsensor[0:MT], ptzsensor[0:MT])
  /* Update Equations */
  for (Time = 0; Time < MaxTime; Time++) {

    CompCurlE(g); 

    TFSF_E(Time, g); 
 
    updateIntE(g); 
  
    updateH(g); 
    
    CompCurlH(g); 

    TFSF_H(Time, g);

    updateIntH(g);

    updateD(g);

    updateE(g, Time);

  } // end of time-stepping

  printf("Done time stepping \n");
   /* Write Sensor data to text files */
  Transmission(g);
  IncSensor(g);
  RefSensor(g);
return 0;
}