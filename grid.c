#include "fdtd-macro.h"
#include "fdtd-alloc.h"
#include "p-macro.h"
#include <complex.h>
#include <math.h>
#include <openacc.h>
// #include <cuda_runtime_api.h>
void gridInit(Grid *g) {
  double imp0 = 377.0, ddx;
  
  int mm, nn, pp, i;
  double XCenter1, XCenter2, YCenter1, YCenter2, r2, XLocC, YLocC, rad, XLoc1, XLoc2, YLoc1, YLoc2, XLoc3, YLoc3, XLoc4, YLoc4, XLoc5, YLoc5, dist1, dist2, dist3, dist4, dist5;
  double epsr, a, PML_no, eps0, c0, dt;
  Type = threeDGrid;   /*@ \label{grid3dhomoA} @*/
  a = 50;

  /* Defining structural constants from p-macro.h */
  PML_no = PMLs;
  MaxTime = MT;   
  SizeX = SizeX0; 
  SizeY = SizeY0;
  SizeZ = SizeZ0;


  c0 = 3e8;
  Cdtds = 1.0/2; // Courant number /*@ \label{grid3dhomoB} @*/
  ddx = DDz0;
  dt = Cdtds*ddx/c0;
  eps0 = 8.85418782e-12;
  
  /* Bloch boundary constants */
  Phix = cexp(I*0*SizeX*ddx);
  Phiy = cexp(I*0*SizeY*ddx);
  Phiz = cexp(I*0*SizeY*ddx);

  printf("SizeX is %d \t SizeY is %d \t SizeZ is %d, dt is %g \n", SizeX0, SizeY0, SizeZ0, dt);
  /* memory allocationi*/
  ALLOC_3D(g->hx,   SizeX, SizeY, SizeZ, complex double); 
  ALLOC_3D(g->hy,   SizeX, SizeY, SizeZ, complex double);
  ALLOC_3D(g->hz,   SizeX, SizeY, SizeZ, complex double);

  ALLOC_3D(g->ex,   SizeX, SizeY, SizeZ, complex double);
  ALLOC_3D(g->ey,   SizeX, SizeY, SizeZ, complex double);
  ALLOC_3D(g->ez,   SizeX, SizeY, SizeZ, complex double);

  ALLOC_3D(g->e1x,  SizeX, SizeY, SizeZ, complex double);
  ALLOC_3D(g->e1y,  SizeX, SizeY, SizeZ, complex double);
  ALLOC_3D(g->e1z,  SizeX, SizeY, SizeZ, complex double);

  ALLOC_3D(g->dx,   SizeX, SizeY, SizeZ, complex double);
  ALLOC_3D(g->dy,   SizeX, SizeY, SizeZ, complex double);
  ALLOC_3D(g->dz,   SizeX, SizeY, SizeZ, complex double);

  /*Curl arrays*/
  ALLOC_3D(g->Curlhx,   SizeX, SizeY, SizeZ, complex double); 
  ALLOC_3D(g->Curlhy,   SizeX, SizeY, SizeZ, complex double); 
  ALLOC_3D(g->Curlhz,   SizeX, SizeY, SizeZ, complex double); 

  ALLOC_3D(g->Curlex,   SizeX, SizeY, SizeZ, complex double); 
  ALLOC_3D(g->Curley,   SizeX, SizeY, SizeZ, complex double); 
  ALLOC_3D(g->Curlez,   SizeX, SizeY, SizeZ, complex double); 

  /* PML memory allocationi*/
  ALLOC_1D(g->sigmaz,    SizeZ, complex double);

  ALLOC_1D(g->PMLhx_0,   SizeZ, complex double);
  ALLOC_1D(g->PMLhx_1,   SizeZ, complex double);
  ALLOC_1D(g->PMLhx_2,   SizeZ, complex double);

  ALLOC_1D(g->PMLhz_2,   SizeZ, complex double);
  ALLOC_1D(g->PMLhz_3,   SizeZ, complex double); 

  ALLOC_1D(g->PMLdx_0,   SizeZ, complex double);
  ALLOC_1D(g->PMLdx_1,   SizeZ, complex double);
  ALLOC_1D(g->PMLdx_2,   SizeZ, complex double);
 
  ALLOC_1D(g->PMLdz_2,   SizeZ, complex double);
  ALLOC_1D(g->PMLdz_3,   SizeZ, complex double);

  /*PML integration arrays*/
  ALLOC_3D(g->IChz,   SizeX, SizeY, SizeZ, complex double); 
  ALLOC_3D(g->ICez,   SizeX, SizeY, SizeZ, complex double); 

  /*memory allocation for z transform fields for X comp*/  

  // ALLOC_3D(g->ixx,   SizeX, SizeY, SizeZ, complex double);

  ALLOC_3D(g->sxx,   SizeX, SizeY, SizeZ, complex double);
  ALLOC_3D(g->s1xx,  SizeX, SizeY, SizeZ, complex double);
  ALLOC_3D(g->s2xx,  SizeX, SizeY, SizeZ, complex double);

  /*memory allocation for z transform fields for Y comp*/
  // ALLOC_3D(g->iyy,   SizeX, SizeY, SizeZ, complex double);

  ALLOC_3D(g->syy,   SizeX, SizeY, SizeZ, complex double);
  ALLOC_3D(g->s1yy,  SizeX, SizeY, SizeZ, complex double);
  ALLOC_3D(g->s2yy,  SizeX, SizeY, SizeZ, complex double);

  /*memory allocation for z transform fields for Z comp*/
  // ALLOC_3D(g->iz,   SizeX, SizeY, SizeZ, complex double);
  ALLOC_3D(g->sz,   SizeX, SizeY, SizeZ, complex double);
  ALLOC_3D(g->s1z,  SizeX, SizeY, SizeZ, complex double);
  ALLOC_3D(g->s2z,  SizeX, SizeY, SizeZ, complex double);
  
  ALLOC_3D(g->epsR, SizeX, SizeY, SizeZ,  double);
  ALLOC_3D(g->media, SizeX, SizeY, SizeZ,  int);

  /* Sensor array memory allocation */
  ALLOC_1D(g->rxsensor,   MaxTime, complex double);
  ALLOC_1D(g->ixsensor,   MaxTime, complex double);
  ALLOC_1D(g->txsensor,   MaxTime, complex double);

  ALLOC_1D(g->rysensor,   MaxTime, complex double);
  ALLOC_1D(g->iysensor,   MaxTime, complex double);
  ALLOC_1D(g->tysensor,   MaxTime, complex double);

  ALLOC_1D(g->rzsensor,   MaxTime, complex double);
  ALLOC_1D(g->izsensor,   MaxTime, complex double);
  ALLOC_1D(g->tzsensor,   MaxTime, complex double);
  
  // ALLOC_3D(g->omega_p, SizeX, SizeY, SizeZ,  double);
  // ALLOC_3D(g->omega_b, SizeX, SizeY, SizeZ,  double);
  // ALLOC_3D(g->nu_c,    SizeX, SizeY, SizeZ,  double);

 
  rad = 17.274; 
  XCenter1 = 0;
  YCenter1 = 0;
  XCenter2 = 65;
  YCenter2 = 37;
  pp = 1;
  r2 = rad*rad;



  /*Media = 1 is magnetized plasma
    Media = 2 is non dispersive
    Media = 3 is reversed biased magnetized plasma*/
  /* Arrays are zeroed when memory is allocated to them*/
  /* Arraays that need to be initialized to other values are reinitialized here*/
  for (mm = 0; mm < SizeX; mm++)
    for (nn = 0; nn < SizeY; nn++) 
       for (pp = 0; pp < SizeZ; pp++) {
          EpsR(mm, nn, pp) = 1.0;
          Media(mm, nn, pp) = 2;
	}


  /* PML conductivity definition only done for 1 dimension as only applied in the Z direction*/
  /* SF PML Sigma*/
  double arg;
  int pmpm;
  for (pp = 0; pp < PML_no + 1; pp++) {
    arg = (pp)/(PML_no);
    pmpm = PML_no - pp;
    SigmaZ(pmpm) = (eps0/(2*dt))*arg*arg*arg;
      // printf("pp is %d, g is %f\n",pp, (PML_no-pp)/(PML_no) );
  }
  /* TF PML Sigma*/
  for (pp = 0; pp < PML_no + 1; pp++) {
    pmpm = (SizeZ - PML_no) + pp;
    arg = (pp)/PML_no;
    SigmaZ(pmpm) = (eps0/(2*dt))*arg*arg*arg;

      // printf("pp is %d, g is %g \n",pp, (1- ((SizeZ - 1 - pp)/(PML_no))));
  }

  /*Dx PML update constants*/
  for (pp = 0; pp < SizeZ; pp++) {
    PMLDx_0(pp) = (1/dt) + (SigmaZ(pp)/(2*eps0));
    PMLDx_1(pp) = ((1/dt) - (SigmaZ(pp)/(2*eps0)))/PMLDx_0(pp);
    PMLDx_2(pp) = (c0/PMLDx_0(pp));
    
  }
  /*Dy PML update constants*/
  /*Not defined as identical Dx terms are reused*/
  // for (mm = 0; mm < SizeX; mm++) 
  //   for (nn = 0; nn < SizeY; nn++) 
  //     for (pp = 0; pp < SizeZ; pp++) {
  //       PMLDy_0(mm,nn,pp) = (1/dt) + (SigmaZ(mm,nn, pp)/(2*eps0));
  //       PMLDy_1(mm,nn,pp) = ((1/dt) - (SigmaZ(mm,nn, pp)/(2*eps0)))/PMLDy_0(mm,nn,pp);
  //       PMLDy_2(mm,nn,pp) = (c0/PMLDy_0(mm,nn,pp));
  //     }
  /*Dz PML update constants*/ 
  for (pp = 0; pp < SizeZ; pp++) {
    PMLDz_2(pp) = c0*dt;
    PMLDz_3(pp) = c0*dt*dt*SigmaZ(pp)/eps0;
  }


  for (pp = 0; pp < SizeZ - 1; pp++) {
    SigmaZ(pp) = (SigmaZ(pp) + SigmaZ(pp + 1))/2;
  }
        
  
  /* PML constants for the Z boundary only PML*/
  /*Hx PML update constants*/
 
  for (pp = 0; pp < SizeZ; pp++) {
    PMLHx_0(pp) =  (1/dt) + (SigmaZ(pp)/(2*eps0));

    PMLHx_1(pp) = ((1/dt) - (SigmaZ(pp)/(2*eps0)))/PMLHx_0(pp);
    printf("PMLHx_1 is %g \n", PMLHx_1(pp));
    /*assuming relative permiability of 1*/
    PMLHx_2(pp) = -1*(c0/PMLHx_0(pp));
    printf("pp is %d, SigmaZ is %g , PMLHx_0 is %g \n", pp, SigmaZ(pp), PMLHx_0(pp));
    printf("PMLHx_2 is %g \n", PMLHx_2(pp));
  }
  /*Hy PML update constants*/
  /*Not defined as identical Hx terms are reused*/ 
  // for (mm = 0; mm < SizeX; mm++) 
  //   for (nn = 0; nn < SizeY; nn++) 
  //     for (pp = 0; pp < SizeZ; pp++) {
  //       PMLHy_0(mm,nn,pp) = (1/dt) + (SigmaZ(mm,nn, pp)/(2*eps0));
  //       PMLHy_1(mm,nn,pp) = ((1/dt) - (SigmaZ(mm,nn, pp)/(2*eps0)))/PMLHy_0(mm,nn,pp);
  //       /*assuming relative perm of 1*/
  //       PMLHy_2(mm,nn,pp) = -1*(c0/PMLHy_0(mm,nn,pp));
  //     }

  /*Hz PML update constants*/ 
  for (pp = 0; pp < SizeZ; pp++) {
    /*assuming relative perm of 1*/
    PMLHz_2(pp) = -1*c0*dt;
    PMLHz_3(pp) = -1*c0*dt*dt*SigmaZ(pp)/eps0;
  }
printf("Hz PML assigned \n");




/* putting in the Ag layer*/
  for (mm = 0; mm < SizeX; mm++)
      for (nn = 0; nn < SizeY; nn++) 
        for (pp = (PML_no + 500); pp <(PML_no + 550); pp++) {
          Media(mm, nn, pp)  = 3;

        }

  /* Hexagonal cylyndrical hole assignment*/
  // for (i = 0; i<31; i++){
  //   for (mm = 0; mm < SizeX; mm++) /*@ \label{grid3dhomoE} @*/
  //     for (nn = 0; nn < SizeY; nn++)
	 //     for (pp = 0; pp < SizeZ; pp++) {

  //     	  XLoc1 = XCenter1 + i*a*sqrt(3) - mm;
  //     	  XLoc2 = XCenter1 + i*a*sqrt(3) - mm;
  //     	  XLoc3 = XCenter1 + i*a*sqrt(3) - mm;
  //     	  XLoc4 = XCenter1 + a*sqrt(3)/2 + i*a*sqrt(3) - mm;
  //     	  XLoc5 = XCenter1 + a*sqrt(3)/2 + i*a*sqrt(3) - mm;

  //     	  YLoc1 = YCenter1 - nn;
  //     	  YLoc2 = YCenter1 + a - nn;
  //     	  YLoc3 = YCenter1 + 2*a - nn;
  //     	  YLoc4 = YCenter1 + a/2 - nn;
  //     	  YLoc5 = YCenter1 + 3*a/2 - nn;

  //     	  dist1 = (int)(sqrt(XLoc1*XLoc1 + YLoc1*YLoc1));
  //     	  dist2 = (int)(sqrt(XLoc2*XLoc2 + YLoc2*YLoc2));
  //     	  dist3 = (int)(sqrt(XLoc3*XLoc3 + YLoc3*YLoc3));
  //     	  dist4 = (int)(sqrt(XLoc4*XLoc4 + YLoc4*YLoc4));
  //     	  dist5 = (int)(sqrt(XLoc5*XLoc5 + YLoc5*YLoc5));
  //     	  if((dist1 < rad)||(dist2 < rad)||(dist3 < rad)||(dist4 < rad)||(dist5 < rad)){
  //           Media(mm, nn, pp) = 2;
  //     	    EpsR(mm, nn, pp) = 14.0;
  //     	  }
  //     	}
  // }


/*media visualization */
  float temp;
  char filename[100];
  FILE *out;
  sprintf(filename,"Media/Media.dat");
  out = fopen(filename,"w");
  mm = SizeX/2;
  //nn = 5;
  for (pp = 0; pp < SizeZ; pp++) {
    for (nn = 0; nn < SizeY; nn++){
      temp = (float)Media(mm,nn, pp); // store data as a float
      // if(pp > PML_no+355){
      //   temp = 4;}
      if(pp < PML_no){
          temp = 0;}
      if(pp> (SizeZ - PML_no)){
          temp = 0;
      }
      
      fprintf(out, "%d \t %f \n", pp, temp);
    }
  }

  fclose(out);

  sprintf(filename,"Media/XY.dat");
  out = fopen(filename,"w");
  pp = PML_no;
  for (nn = 0; nn < SizeY; nn++) {
    for (mm=0; mm<SizeX; mm++) {
      temp = (float)Media(mm,nn, pp); // store data as a float
      fprintf(out, "%f \n", temp);
    }
  }
  fclose(out);

  return;
}  /* end gridInit() */
