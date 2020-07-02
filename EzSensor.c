#include <stdio.h>
#include <stdlib.h>
#include "fdtd-macro.h"
#include "p-macro.h"
#include <complex.h>
#include <math.h>
#include <openacc.h>

/*code to record a sensor of the data*/
void sensorInit(Grid *g){
	char filename[100];
	FILE *out;
	sprintf(filename, "Trans/Trans.txt");
        out = fopen(filename, "w");
	fprintf(out, "Trans init \n");
	fclose(out);

        sprintf(filename, "Inc/Inc.txt");
        out = fopen(filename, "w");
        fprintf(out, "Inc Sensor\n");
        fclose(out);

        sprintf(filename, "Ref/Ref.txt");
        out = fopen(filename, "w");
        fprintf(out, "Ref Sensor \n");
        fclose(out);


return;
}

void IncSensor(Grid *g){
        int t;
        complex double *pixsensor = g->ixsensor;
        complex double *piysensor = g->iysensor;
        complex double *pizsensor = g->izsensor;
        #pragma acc exit data copyout(pixsensor[0:MT], piysensor[0:MT], pizsensor[0:MT])
        double tempX, tempY, tempZ, itempX, itempY, itempZ;
        char filename[100];
        FILE *out;
        sprintf(filename, "Inc/Inc.txt");
        out = fopen(filename, "a");
        for(t = 0; t < MT; t++){              
                /* print the time stamp and the Ex field right before the QWS*/
                tempX = creal(IXSensor(t));
                tempY = creal(IYSensor(t));
                tempZ = creal(IZSensor(t));

                itempX = cimag(IXSensor(t));
                itempY = cimag(IYSensor(t));
                itempZ = cimag(IZSensor(t));
                fprintf(out, "%g \t %g \t %g \t %g \t %g \t %g \n", tempX, itempX, tempY, itempY, tempZ, itempZ);//, cimag(Hz(80, 34, 1)), creal(Hz(33, 43, 1)), cimag(Hz(33, 43, 1)), creal(Hz(83, 8, 1)), cimag(Hz(83, 8, 1)), creal(Hz(81, 44, 1)), cimag(Hz(81, 44, 1)), creal(Hz(36, 44, 1)), cimag(Hz(36, 44, 1)));
        }
        fclose(out);

}


void RefSensor(Grid *g){
        int t;
        complex double *prxsensor = g->rxsensor;
        complex double *prysensor = g->rysensor;
        complex double *przsensor = g->rzsensor;
        #pragma acc exit data copyout(prxsensor[0:MT], prysensor[0:MT], przsensor[0:MT])
        double tempX, tempY, tempZ, itempX, itempY, itempZ;
        char filename[100];
        FILE *out;
        sprintf(filename, "Ref/Ref.txt");
        out = fopen(filename, "a");
        for(t = 0; t < MT; t++){              
                /* print the time stamp and the Ex field right before the QWS*/
                tempX = creal(RXSensor(t));
                tempY = creal(RYSensor(t));
                tempZ = creal(RZSensor(t));

                itempX = cimag(RXSensor(t));
                itempY = cimag(RYSensor(t));
                itempZ = cimag(RZSensor(t));
                fprintf(out, "%g \t %g \t %g \t %g \t %g \t %g \n", tempX, itempX, tempY, itempY, tempZ, itempZ);//, cimag(Hz(80, 34, 1)), creal(Hz(33, 43, 1)), cimag(Hz(33, 43, 1)), creal(Hz(83, 8, 1)), cimag(Hz(83, 8, 1)), creal(Hz(81, 44, 1)), cimag(Hz(81, 44, 1)), creal(Hz(36, 44, 1)), cimag(Hz(36, 44, 1)));
        }
        fclose(out);

return;      
}

void Transmission(Grid *g){
        int t;
        complex double *ptxsensor = g->txsensor;
        complex double *ptysensor = g->tysensor;
        complex double *ptzsensor = g->tzsensor;
        #pragma acc exit data copyout(ptxsensor[0:MT], ptysensor[0:MT], ptzsensor[0:MT])
        double tempX, tempY, tempZ, itempX, itempY, itempZ;
        char filename[100];
        FILE *out;
        sprintf(filename, "Trans/Trans.txt");
        out = fopen(filename, "a");
        for(t = 0; t < MT; t++){              
                /* print the time stamp and the Ex field right before the QWS*/
                tempX = creal(TXSensor(t));
                tempY = creal(TYSensor(t));
                tempZ = creal(TZSensor(t));

                itempX = cimag(TXSensor(t));
                itempY = cimag(TYSensor(t));
                itempZ = cimag(TZSensor(t));
                fprintf(out, "%g \t %g \t %g \t %g \t %g \t %g \n", tempX, itempX, tempY, itempY, tempZ, itempZ);//, cimag(Hz(80, 34, 1)), creal(Hz(33, 43, 1)), cimag(Hz(33, 43, 1)), creal(Hz(83, 8, 1)), cimag(Hz(83, 8, 1)), creal(Hz(81, 44, 1)), cimag(Hz(81, 44, 1)), creal(Hz(36, 44, 1)), cimag(Hz(36, 44, 1)));
        }
        fclose(out);

return;      
}




