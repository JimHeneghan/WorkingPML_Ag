#include "fdtd-macro.h"
#include "p-macro.h"
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <openacc.h>




#define SQR(x) ((x)*(x))

/* update magnetic field */

void CompCurlE(Grid *g){
  double DDx = DDx0;
  double DDy = DDy0;
  double DDz = DDz0;
  // printf("started CurlE \n");
  /*Calculate the Hx field*/
  int mm, nn, pp;

  int tc = SizeX0*SizeY0*SizeZ0;
  complex double *pCurlex = g->Curlex;
  complex double *pCurley = g->Curley;
  complex double *pCurlez = g->Curlez;

  complex double *pex = g->ex;
  complex double *pey = g->ey;
  complex double *pez = g->ez;
  // #pragma acc parallel  
  // {
        // #pragma acc loop independent collapse(3)
  #pragma acc kernels loop independent collapse(3) present(pCurlex[0:tc], pey[0:tc],pez[0:tc]) 
    for (mm = 0; mm < SizeX0; mm++){
      // #pragma acc loop independent
      for (nn = 0; nn < SizeY0 - 1; nn++){
        // #pragma acc loop independent
        for (pp = 0; pp < SizeZ0 - 1; pp++){
          pCurlEx(mm, nn, pp) = (pEz(mm, nn + 1, pp) - pEz(mm, nn, pp))/DDy
            - (pEy(mm, nn, pp + 1) - pEy(mm, nn, pp))/DDz;
      }
     }
    }

  // #pragma acc kernels loop independent collapse(2) present(pCurlex[0:tc], pey[0:tc],pez[0:tc])  
  //   for (mm = 0; mm < SizeX0; mm++){
  //     for (nn = 0; nn < SizeY0 - 1; nn++){
  //       pCurlEx(mm, nn, SizeZ0-1) = (pEz(mm, nn + 1, SizeZ0-1) - pEz(mm, nn, SizeZ0-1))/DDy
  //         - (0 - pEy(mm, nn, SizeZ0-1))/DDz;
  //    }
  //   }

  #pragma acc kernels loop independent collapse(2) present(pCurlex[0:tc], pey[0:tc],pez[0:tc])  
    for (mm = 0; mm < SizeX0; mm++){
      for (pp = 0; pp < SizeZ0 -1; pp++){
       pCurlEx(mm, SizeY0 - 1, pp) = (pEz(mm, 0, pp) - pEz(mm, SizeY0 - 1, pp))/DDy
        - (pEy(mm, SizeY0 - 1, pp + 1) - pEy(mm, SizeY0 -1, pp))/DDz;
      }
    }

  // #pragma acc kernels loop independent present(pCurlex[0:tc], pey[0:tc],pez[0:tc])  
  //   for (mm = 0; mm < SizeX0; mm++){
  //      pCurlEx(mm, SizeY0 - 1, SizeZ0 - 1) = (pEz(mm, 0, SizeZ0 - 1) - pEz(mm, SizeY0 - 1, SizeZ0 - 1))/DDy
  //       - (0 - pEy(mm, SizeY0 -1, SizeZ0 - 1))/DDz;
  //   }

    /* Calculate the Curl Ey field */
  #pragma acc kernels loop independent collapse(3) present(pCurley[0:tc], pex[0:tc],pez[0:tc]) 
    for (mm = 0; mm < SizeX0 - 1; mm++){
      for (nn = 0; nn < SizeY0; nn++){
        for (pp = 0; pp < SizeZ0 -1; pp++){
         pCurlEy(mm, nn, pp) = (pEx(mm, nn, pp + 1) - pEx(mm, nn, pp))/DDz
          - (pEz(mm + 1, nn, pp) - pEz(mm, nn, pp))/DDx;
        }
      }
   }

  // #pragma acc kernels loop independent collapse(2) present(pCurley[0:tc], pex[0:tc],pez[0:tc]) 
  //  for (mm = 0; mm < SizeX0 - 1; mm++){
  //     for (nn = 0; nn < SizeY0; nn++){
  //        pCurlEy(mm, nn, SizeZ0-1) = (0 - pEx(mm, nn, SizeZ0-1))/DDz
  //         - (pEz(mm + 1, nn, SizeZ0-1) - pEz(mm, nn, SizeZ0-1))/DDx;
  //     }
  //  }

  #pragma acc kernels loop independent collapse(2) present(pCurley[0:tc], pex[0:tc],pez[0:tc]) 
    for (nn = 0; nn < SizeY0; nn++){
       for (pp = 0; pp < SizeZ0 -1; pp++){
         pCurlEy(SizeX0 - 1, nn, pp) = (pEx(SizeX0 - 1, nn, pp + 1) - pEx(SizeX0 - 1, nn, pp))/DDz
          - (pEz(0, nn, pp) - pEz(SizeX0 - 1, nn, pp))/DDx;
        }
      }

  // #pragma acc kernels loop independent present(pCurley[0:tc], pex[0:tc],pez[0:tc])   
  //   for (nn = 0; nn < SizeY0; nn++){
  //     pCurlEy(SizeX0 - 1, nn, SizeZ0-1) = (0 - pEx(SizeX0 - 1, nn, SizeZ0-1))/DDz
  //       - (pEz(0, nn, SizeZ0-1) - pEz(SizeX0 - 1, nn, SizeZ0-1))/DDx;
  //   }

      /* Calculate the Curl of Ez field */ 
  #pragma acc kernels loop independent collapse(3) present(pCurlez[0:tc], pex[0:tc],pey[0:tc])     
    for (mm = 0; mm < SizeX0 - 1; mm++){
     for (nn = 0; nn < SizeY0 - 1; nn++){
       for (pp = 0; pp < SizeZ0 -1; pp++){
         pCurlEz(mm, nn, pp) = (pEy(mm + 1, nn, pp) - pEy(mm, nn, pp))/DDx
          - (pEx(mm, nn + 1, pp) - pEx(mm, nn, pp))/DDy;
        }
      }
    }
    #pragma acc kernels loop independent collapse(2) present(pCurlez[0:tc], pex[0:tc],pey[0:tc])  
     for (nn = 0; nn < SizeY0 - 1; nn++){
       for (pp = 0; pp < SizeZ0 -1; pp++){
         pCurlEz(SizeX0 - 1, nn, pp) = (pEy(0, nn, pp) - pEy(SizeX0 - 1, nn, pp))/DDx
          - (pEx(SizeX0 - 1, nn + 1, pp) - pEx(SizeX0 - 1, nn, pp))/DDy; 
        }
      }
  #pragma acc kernels loop independent collapse(2) present(pCurlez[0:tc], pex[0:tc],pey[0:tc]) 
   for (mm = 0; mm < SizeX0 - 1; mm++){
     for (pp = 0; pp < SizeZ0 -1; pp++){
       pCurlEz(mm, SizeY0  - 1, pp) = (pEy(mm + 1, SizeY0  - 1, pp) - pEy(mm, SizeY0  - 1, pp))/DDx
        - (pEx(mm, 0, pp) - pEx(mm, SizeY0  - 1, pp))/DDy;
     }
   }
  #pragma acc kernels loop independent present(pCurlez[0:tc], pex[0:tc],pey[0:tc]) 
   for (pp = 0; pp < SizeZ0 -1; pp++){
     pCurlEz(SizeX0 - 1, SizeY0 - 1, pp) = (pEy(0, SizeY0 - 1, pp) - pEy(SizeX0 - 1, SizeY0 - 1, pp))/DDx
      - (pEx(SizeX0 - 1, 0, pp) - pEx(SizeX0 - 1, SizeY0  - 1, pp))/DDy;
   }
// }
// #pragma acc exit data copyout(pex[0:tc], pey[0:tc], pez[0:tc], pCurlex[0:tc], pCurley[0:tc], pCurlez[0:tc]) 
   // }

  #pragma acc wait
  return;
} /*End Comp Curl E*/



void updateIntE(Grid *g){
  int mm, nn, pp;
  int tc = SizeX0*SizeY0*SizeZ0;
  complex double *pCurlez = g->Curlez;
  complex double *pCurlex = g->Curlex;
  complex double *pCurley = g->Curley;
  complex double *pICez = g->ICez;
  // #pragma acc enter data copyin(pCurlex[0:tc], pCurley[0:tc], pCurlez[0:tc]) 
  {
    #pragma acc kernels loop independent collapse(3) present(pICez[:tc], pCurlez[0:tc]) 
    for (mm = 0; mm < SizeX0; mm++){
       for (nn = 0; nn < SizeY0; nn++){
         for (pp = 0; pp < SizeZ0 - 1; pp++){
          pICEz(mm, nn, pp) = pICEz(mm, nn, pp) + pCurlEz(mm, nn, pp);
         }
       }
    }
  }
   return;
}



void updateH(Grid *g) {
  int mm, nn, pp;
  int tc = SizeX0*SizeY0*SizeZ0;
  /*create new pointers that point to the struct pointers*/
  /*create new macro's for those pointers*/
  /*profit*/
  complex double *pCurlex = g->Curlex;
  complex double *pCurley = g->Curley;
  complex double *pCurlez = g->Curlez;
  complex double *pICez = g->ICez;

  complex double *pPMLhx_1 = g->PMLhx_1;
  complex double *pPMLhx_2 = g->PMLhx_2;

  complex double *pPMLhz_2 = g->PMLhz_2;
  complex double *pPMLhz_3 = g->PMLhz_3;
  complex double *phx = g->hx;
  complex double *phy = g->hy;
  complex double *phz = g->hz;
    /*Calculate the Hx field*/
    #pragma acc parallel loop async independent collapse(3) present(pCurlex[0:tc], phx[0:tc], pPMLhx_1[0:SizeZ0], pPMLhx_2[0:SizeZ0])
    for (mm = 0; mm < SizeX0; mm++){
      for (nn = 0; nn < SizeY0; nn++){
        for (pp = 0; pp < SizeZ0 - 1; pp++){
          pHx(mm, nn, pp) =  pHx(mm, nn, pp) * pPMLHx_1(pp)
            + pCurlEx(mm, nn, pp)*pPMLHx_2(pp);
        }
      }
    }

    // #pragma acc exit data copyout(pCurlex[0:tc], pPMLhx_2[0:SizeZ0])
    // // for (pp = 0; pp < SizeZ0 - 1; pp++){
    //   printf("CurlEx is %g\n",CurlEx(5, 5, 130));
    // // }
    // #pragma acc enter data copyin(pCurlex[0:tc], pPMLhx_2[0:SizeZ0])

    /* Calculate the Hy field */    
  #pragma acc parallel loop async independent collapse(3) present(pCurley[0:tc],phy[0:tc],pPMLhx_1[0:SizeZ0],pPMLhx_2[0:SizeZ0]) 
    for (mm = 0; mm < SizeX0; mm++){
      for (nn = 0; nn < SizeY0; nn++){
        for (pp = 0; pp < SizeZ0 - 1; pp++){
         pHy(mm, nn, pp) =  pHy(mm, nn, pp) * pPMLHx_1(pp) 
          + pCurlEy(mm, nn, pp) * pPMLHx_2(pp);
        }
      }
    }

   /* Calculate the Hz field */
  #pragma acc parallel loop async independent collapse(3) present(pCurlez[0:tc],phz[0:tc],pPMLhz_3[0:SizeZ0],pPMLhz_2[0:SizeZ0],pICez[:tc])     
   for (mm = 0; mm < SizeX0; mm++){
     for (nn = 0; nn < SizeY0; nn++){
       for (pp = 0; pp < SizeZ0 - 1; pp++){
         pHz(mm, nn, pp) = pHz(mm, nn, pp) 
          + pCurlEz(mm, nn, pp) * pPMLHz_2(pp)
            + pICEz(mm, nn, pp) * pPMLHz_3(pp); 
       }
     }
   }

  #pragma acc wait
  return;
}  /* end updateH() */


/* Calculate Curl Hx field */
void CompCurlH(Grid *g){
    
  int mm, nn, pp;
  int tc = SizeX0*SizeY0*SizeZ0;

  complex double *pCurlhx = g->Curlhx;
  complex double *pCurlhy = g->Curlhy;
  complex double *pCurlhz = g->Curlhz;

  complex double *phx = g->hx;
  complex double *phy = g->hy;
  complex double *phz = g->hz;
  double DDx = DDx0;
  double DDy = DDy0;
  double DDz = DDz0;

  #pragma acc kernels loop independent collapse(3) present(pCurlhx[0:tc], phy[0:tc],phz[0:tc]) 
    for (mm = 0; mm < SizeX0; mm++){
      for (nn = 1; nn < SizeY0; nn++){
        for (pp = 1; pp < SizeZ0; pp++){
          pCurlHx(mm, nn, pp) = (pHz(mm, nn, pp) - pHz(mm, nn - 1, pp))/DDy -
            (pHy(mm, nn, pp) - pHy(mm, nn, pp - 1))/DDz;
        }
      }
    }

  // #pragma acc kernels loop independent collapse(2) present(pCurlhx[0:tc], phy[0:tc],phz[0:tc])
  //   for (mm = 0; mm < SizeX0; mm++){
  //     for (nn = 1; nn < SizeY0; nn++){
  //       pCurlHx(mm, nn, 0) = (pHz(mm, nn, 0) - pHz(mm, nn - 1, 0))/DDy -
  //         (pHy(mm, nn, 0) - 0)/DDz;
  //     }
  //   }

  #pragma acc kernels loop independent collapse(2) present(pCurlhx[0:tc], phy[0:tc],phz[0:tc]) 
    for (mm = 0; mm < SizeX0; mm++){
      for (pp = 1; pp < SizeZ0; pp++){
        pCurlHx(mm, 0, pp) = (pHz(mm, 0, pp) - pHz(mm, SizeY0 - 1, pp))/DDy -
         (pHy(mm, 0, pp) - pHy(mm, 0, pp - 1))/DDz;
      }
    }

  // #pragma acc kernels loop independent present(pCurlhx[0:tc], phy[0:tc],phz[0:tc])
  //   for (mm = 0; mm < SizeX0; mm++){
  //       pCurlHx(mm, 0, 0) = (pHz(mm, 0, 0) - pHz(mm, SizeY0 - 1, 0))/DDy -
  //         (pHy(mm, 0, 0) - 0)/DDz;      
  //   }

  /* Calculate Curl Hy field */
  #pragma acc kernels loop independent collapse(3) present(pCurlhy[0:tc], phx[0:tc],phz[0:tc])     
    for (mm = 1; mm < SizeX0; mm++){
      for (nn = 0; nn < SizeY0; nn++){
        for (pp = 1; pp < SizeZ0; pp++){
          pCurlHy(mm, nn, pp) = (pHx(mm, nn, pp) - pHx(mm, nn, pp - 1))/DDz -
            (pHz(mm, nn, pp) - pHz(mm - 1, nn, pp))/DDx;
        }
      }
    }

  // #pragma acc kernels loop independent collapse(2) present(pCurlhy[0:tc], phx[0:tc],phz[0:tc])     
  //   for (nn = 0; nn < SizeY0; nn++){
  //     for (pp = 1; pp < SizeZ0; pp++){
  //       pCurlHy(mm, nn, 0) = (pHx(mm, nn, 0) - 0)/DDz -
  //         (pHz(mm, nn, 0) - pHz(mm - 1, nn, 0))/DDx;
  //     }
  //   }

  #pragma acc kernels loop independent collapse(2) present(pCurlhy[0:tc], phx[0:tc],phz[0:tc]) 
    for (nn = 0; nn < SizeY0; nn++){
      for (pp = 1; pp < SizeZ0; pp++){
        pCurlHy(0, nn, pp) = (pHx(0, nn, pp) - pHx(0, nn, pp - 1))/DDz -
          (pHz(0, nn, pp) - pHz(SizeX0 - 1, nn, pp))/DDx;
      }
    }

  // #pragma acc kernels loop independent present(pCurlhy[0:tc], phx[0:tc],phz[0:tc]) 
  //   for (nn = 0; nn < SizeY0; nn++){
  //     pCurlHy(0, nn, 0) = (pHx(0, nn, 0) - 0)/DDz -
  //       (pHz(0, nn, 0) - pHz(SizeX0 - 1, nn, 0))/DDx;
  //   }

  /* Calculate Curl Hz field */
  #pragma acc kernels loop independent collapse(3) present(pCurlhz[0:tc], phx[0:tc],phy[0:tc]) 
    for (mm = 1; mm < SizeX0; mm++){
      for (nn = 1; nn < SizeY0; nn++){
        for (pp = 1; pp < SizeZ0; pp++){
          pCurlHz(mm, nn, pp) = (pHy(mm, nn, pp) - pHy(mm - 1, nn, pp))/DDx -
            (pHx(mm, nn, pp) - pHx(mm, nn - 1, pp))/DDy;
        }
      }
    }

  #pragma acc kernels loop independent collapse(2) present(pCurlhz[0:tc], phx[0:tc],phy[0:tc]) 
    for (nn = 1; nn < SizeY0; nn++){
      for (pp = 1; pp < SizeZ0; pp++){
        pCurlHz(0, nn, pp) = (pHy(0, nn, pp) - pHy(SizeX0 - 1, nn, pp))/DDx -
          (pHx(0, nn, pp) - pHx(0, nn - 1, pp))/DDy;
      }
    }

  #pragma acc kernels loop independent collapse(2) present(pCurlhz[0:tc], phx[0:tc],phy[0:tc]) 
    for (mm = 1; mm < SizeX0; mm++){
      for (pp = 1; pp < SizeZ0; pp++){ 
        pCurlHz(mm, 0, pp) = (pHy(mm, 0, pp) - pHy(mm - 1, 0, pp))/DDx -
          (pHx(mm, 0, pp) - pHx(mm, SizeY0 - 1, pp))/DDy;
      }
    }

  #pragma acc kernels loop independent present(pCurlhz[0:tc], phx[0:tc],phy[0:tc]) 
    for (pp = 1; pp < SizeZ0; pp++){
       pCurlHz(0, 0, pp) = (pHy(0, 0, pp) - pHy(SizeX0 - 1, 0, pp))/DDx - 
        (pHx(0, 0, pp) - pHx(0, SizeY0 - 1, pp))/DDy;
    }

  #pragma acc wait
  return;
} /*end CompCurlH()*/


/*update H integrations*/
void updateIntH(Grid *g){
  int mm, nn, pp;
  int tc = SizeX0*SizeY0*SizeZ0;

  complex double *pCurlhx = g->Curlhx;
  complex double *pCurlhy = g->Curlhy;
  complex double *pCurlhz = g->Curlhz;

  complex double *pIChz = g->IChz;

    #pragma acc kernels loop independent collapse(3) present(pIChz[:tc], pCurlhz[0:tc]) 
    for (mm = 0; mm < SizeX0; mm++){
      for (nn =0; nn < SizeY0; nn++){
        for (pp = 1; pp < SizeZ0; pp++){
          pICHz(mm, nn, pp) = pICHz(mm, nn, pp) + pCurlHz(mm, nn, pp);
        }
      }
    }

  return;
}



/* update D field */
void updateD(Grid *g) {
  int mm, nn, pp;
  complex double CurlH;
  complex double Psix = 1;//conj(Phix);
  complex double Psiy = 1;//conj(Phiy);
  complex double Psiz = 1;//conj(Phiz);

  int tc = SizeX0*SizeY0*SizeZ0;

  complex double *pCurlhx = g->Curlhx;
  complex double *pCurlhy = g->Curlhy;
  complex double *pCurlhz = g->Curlhz;

  complex double *pdx = g->dx;
  complex double *pdy = g->dy;
  complex double *pdz = g->dz;

  complex double *pPMLdx_1 = g->PMLdx_1;
  complex double *pPMLdx_2 = g->PMLdx_2;

  complex double *pPMLdz_2 = g->PMLdz_2;
  complex double *pPMLdz_3 = g->PMLdz_3;

  complex double *pIChz = g->IChz;

    /* Calculate Dx field */
  #pragma acc parallel loop async independent collapse(3) present(pCurlhx[0:tc], pdx[0:tc], pPMLdx_1[0:SizeZ0], pPMLdx_2[0:SizeZ0]) 
    for (mm = 0; mm < SizeX0; mm++){
      for (nn = 0; nn < SizeY0; nn++){
        for (pp = 1; pp < SizeZ0; pp++){
          pDx(mm, nn, pp) = pDx(mm, nn, pp) * pPMLDx_1(pp) 
             + pCurlHx(mm, nn, pp) * pPMLDx_2(pp);
        }
      }
    }

  /* Calculate Dy field */
  #pragma acc parallel loop async independent collapse(3) present(pCurlhy[0:tc], pdy[0:tc], pPMLdx_1[0:SizeZ0], pPMLdx_2[0:SizeZ0]) 
    for (mm = 0; mm < SizeX0; mm++){
      for (nn = 0; nn < SizeY0; nn++){
        for (pp = 1; pp < SizeZ0; pp++){
          pDy(mm, nn, pp) = pDy(mm, nn, pp) * pPMLDx_1(pp) 
             + pCurlHy(mm, nn, pp) * pPMLDx_2(pp);
        }
      }
    }

  /* Calculate Dz field */
  #pragma acc parallel loop async independent collapse(3) present(pCurlhz[0:tc], pdz[0:tc], pPMLdz_3[0:SizeZ0], pPMLdz_2[0:SizeZ0], pIChz[:tc]) 
    for (mm = 0; mm < SizeX0; mm++){
      for (nn = 0; nn < SizeY0; nn++){
        for (pp = 1; pp < SizeZ0; pp++){
          pDz(mm, nn, pp) = pDz(mm, nn, pp) 
            + pCurlHz(mm, nn, pp) * pPMLDz_2(pp)
              + pPMLDz_3(pp) * pICHz(mm, nn, pp);
        }
      }
    }

  #pragma acc wait

  return;
}  /* end updateD() */


void updateE(Grid *g, int t) {

  int mm, nn, pp;
  /* Declaring Lorentz Equation Constants for hBN*/
  double w_nuXY, S_nu_XY, gam,  eps_inf,   A,  B,  C,  EXP1,  EXP2,  CS,  SN, C1XY, C2XY;
  double w_nuZ,  S_nu_Z,  gamZ, eps_infz,  Az, Bz, Cz, EXP1Z, EXP2Z, CSZ, SNZ, C1Z, C2Z;
  /* Declaring Drude Equation Constants for Ag*/
  double wp, gamAg, EXPAg, eps_inf_Ag, AAg, C1Ag, C2Ag;

  /* Declaring more generic constants*/
  eps_inf = 4.87;
  eps_infz = 2.95;
  double epsilon0 = 8.85418782e-12;
  double dt = DDz0/(2*3e8);
  double Epsr = 1.0;
  /* Redeclaring constants for use in update equations to not be part of a struct */
  int MaxTime0 = MaxTime;
  int PMLs0 = PMLs;
  int t0 = t;
  int tc = SizeX0*SizeY0*SizeZ0;

  /* Declaring array pointers to pointers to Struct:Grid type pointes*/
  complex double *pex = g->ex;
  complex double *pey = g->ey;
  complex double *pez = g->ez;

  complex double *pdx = g->dx; 
  complex double *pdy = g->dy;
  complex double *pdz = g->dz;

  complex double *psxx = g->sxx;
  complex double *psyy = g->syy;
  complex double  *psz = g->sz;

  complex double *ps1xx = g->s1xx;
  complex double *ps1yy = g->s1yy;
  complex double  *ps1z = g->s1z; 

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

  /*Calculate Ex Field*/
  #pragma acc parallel loop async independent collapse(3) present(pex[:tc], pdx[:tc], psxx[:tc], ps1xx[:tc], pmedia[:tc]) 
   for (mm = 0; mm < SizeX0; mm++){
     for (nn = 0; nn < SizeY0; nn++){
       for (pp = 1; pp < SizeZ0; pp++){
        /* If statements to indicate which update equation to use
        /* Media = 1 is Lorentz model hBN
        /* Media = 2 is non dispersive
        /* Media = 3 is Drude Model Ag */
        if (pMedia(mm, nn, pp) == 1){
          pEx(mm, nn, pp)  = (pDx(mm, nn, pp) - pSxx(mm, nn, pp))/eps_inf;
        } else if(pMedia(mm, nn, pp) == 2) {
          pEx(mm, nn, pp)  = pDx(mm, nn, pp)/Epsr;
        } else if(pMedia(mm, nn, pp) == 3) {
          pEx(mm, nn, pp)  = (pDx(mm, nn, pp) - pSxx(mm, nn, pp));
        }
       }
      }
    }
      /*Incident, REflection and Transmission sensors*/
      /*These could be modified to FFt sensors       */
      /*or stand aloneor FFt functions in future     */


  /*Calculate Ey Field*/ 
  #pragma acc parallel loop async independent collapse(3) present(pey[:tc], pdy[:tc], psyy[:tc], ps1yy[:tc], pmedia[:tc])//, prysensor[:MT], piysensor[:MT], ptysensor[:MT]) 
    for (mm = 0; mm < SizeX0; mm++){
      for (nn = 0; nn < SizeY0; nn++){
        for (pp = 1; pp < SizeZ0; pp++){
          if (pMedia(mm, nn, pp) == 1){
            pEy(mm, nn, pp)  = (pDy(mm, nn, pp) -  pSyy(mm, nn, pp))/eps_inf;
          } else if(pMedia(mm, nn, pp) == 2) { 
            pEy(mm, nn, pp)  = pDy(mm, nn, pp)/Epsr;
          } else if(pMedia(mm, nn, pp) == 3) {
            pEy(mm, nn, pp)  = (pDy(mm, nn, pp) - pSyy(mm, nn, pp));// + pS1yy(mm, nn, pp))/eps_inf_Ag;
          }
        }
      }
      // pIYSensor(t) = (pEy(5, 5, 130) + pEy(5, 4, 130))/2;
      // // printf("updated I\n");
      // pRYSensor(t) = (pEy(5, 5, 120) + pEy(5, 4, 120))/2;
      // // printf("updated t\n");
      // pTYSensor(t) = (pEy(5, 5, 245) + pEy(5, 4, 245))/2;
    }

  /*Calculate Ez Field*/
  #pragma acc parallel loop async independent collapse(3) present(pez[:tc], pdz[:tc], psz[:tc], ps1z[:tc], pmedia[:tc])//, przsensor[:MT], pizsensor[:MT], ptzsensor[:MT]) 
    for (mm = 0; mm < SizeX0; mm++){
      for (nn = 0; nn < SizeY0; nn++){
        for (pp = 1; pp < SizeZ0; pp++){
          if (pMedia(mm, nn, pp) == 1){
            pEz(mm, nn, pp) = (pDz(mm, nn, pp) - pSz(mm, nn, pp))/eps_infz;
          } else if(pMedia(mm, nn, pp) == 2) {
            pEz(mm, nn, pp)  = pDz(mm, nn, pp)/Epsr;//(mm, nn, pp);
          } else if(pMedia(mm, nn, pp) == 3) {
            pEz(mm, nn, pp)  = (pDz(mm, nn, pp) - pSz(mm, nn, pp));
        }
      }
      // pIZSensor(t) = pEz(5, 5, 130);
      // // printf("updated I\n");
      // pRZSensor(t) = pEz(5, 5, 120);
      // // printf("updated t\n");
      // pTZSensor(t) = pEz(5, 5, 245);
    }
  }
#pragma acc wait
#pragma acc parallel loop  independent present(prxsensor[:MT], pixsensor[:MT], ptxsensor[:MT], pex[:tc])
  for (pp = 0; pp < 1; pp++){
    pIXSensor(t) = (pEx(5, 5, 20 + PMLs0));
    pRXSensor(t) = (pEx(5, 5, 10 + PMLs0));
    pTXSensor(t) = (pEx(5, 5, 590 + PMLs0));
  }
  /*Lorentz term update equations*/
  /*could be done in the E* update loops but done here to increase portability to an off diagnal medium*/

  /*XY hBN Lorentz constants*/  
    // w_nuXY = 2.58276e14;
    // gam = 6.02526055e12;
    // S_nu_XY = 1.83;

    // A = (gam/2);
    // B = w_nuXY*(sqrt(1 - SQR(gam/(2*w_nuXY))));
    // C = w_nuXY/(sqrt(1 - SQR(gam/(2*w_nuXY))));
    // EXP1 = exp(-A*dt);
    // EXP2 = exp(-2*A*dt);
    // CS = cos(B*dt);
    // SN = sin(B*dt);

    // C1XY = 2*EXP1*CS;
    // C2XY = EXP1*SN*C*dt*S_nu_XY;


  /*Z hBN Lorentz constants*/ 
    // w_nuZ = 1.41e14;
    // gamZ = 3.7891e11;//2.3864e12;
    // S_nu_Z= 0.61;

    // Az = (gamZ/2);
    // Bz = w_nuZ*(sqrt(1 - SQR(gamZ/(2*w_nuZ))));
    // Cz = w_nuZ/(sqrt(1 - SQR(gamZ/(2*w_nuZ))));
    // EXP1Z = exp(-Az*dt);
    // EXP2Z = exp(-2*Az*dt);
    // CSZ = cos(Bz*dt);
    // SNZ = sin(Bz*dt);

    // C1Z = 2*EXP1Z*CSZ;
    // C2Z = EXP1Z*SNZ*Cz*dt*S_nu_Z;

    /*Drude model for Ag from Palik */

    complex double *ps2xx = g->s2xx;
    complex double *ps2yy = g->s2yy;
    complex double  *ps2z = g->s2z;

    complex double *pe1x = g->e1x;
    complex double *pe1y = g->e1y;
    complex double *pe1z = g->e1z;   

    gamAg = 9.79125662e13;
    wp = 1.15136316e16;
    EXPAg =exp(-gamAg*dt);
    AAg = SQR(wp)*dt/gamAg;
    /* Z Transfrom update arrays */
    #pragma acc kernels loop independent collapse(3) present(pex[:tc], pey[:tc], pez[:tc], pe1x[:tc], pe1y[:tc], pe1z[:tc], psxx[:tc], psyy[:tc], psz[:tc], ps1xx[:tc], ps1yy[:tc], ps1z[:tc], ps2xx[:tc], ps2yy[:tc], ps2z[:tc], pmedia[:tc]) 
    for (mm = 0; mm < SizeX0; mm++){
      for (nn = 0; nn < SizeY0; nn++){
        for (pp = 1; pp < SizeZ0; pp++){
          if (pMedia(mm, nn, pp) == 3){

            pSxx(mm, nn, pp) = (1 + EXPAg)*pS1xx(mm, nn, pp) - EXPAg*pS2xx(mm, nn, pp) + AAg*(1 - EXPAg)*pEx(mm, nn, pp);

            pSyy(mm, nn, pp) = (1 + EXPAg)*pS1yy(mm, nn, pp) - EXPAg*pS2yy(mm, nn, pp) + AAg*(1 - EXPAg)*pEy(mm, nn, pp);

            pSz(mm, nn, pp) = (1 + EXPAg)*pS1z(mm, nn, pp) - EXPAg*pS2z(mm, nn, pp) + AAg*(1 - EXPAg)*pEz(mm, nn, pp);

            pS2xx(mm, nn, pp) = pS1xx(mm, nn, pp);
            pS1xx(mm, nn, pp) = pSxx(mm, nn, pp);

            pS2yy(mm, nn, pp) = pS1yy(mm, nn, pp);
            pS1yy(mm, nn, pp) = pSyy(mm, nn, pp);

            pS2z(mm, nn, pp) = pS1z(mm, nn, pp);
            pS1z(mm, nn, pp) = pSz(mm, nn, pp);

          // pE1x(mm, nn, pp) = pEx(mm, nn, pp);
          // pE1y(mm, nn, pp) = pEy(mm, nn, pp);
          // pE1z(mm, nn, pp) = pEz(mm, nn, pp);
          }
        }
      }
    }




return;
}