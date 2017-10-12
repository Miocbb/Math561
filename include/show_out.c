#include<stdio.h>
#include"matrix.h"
#include"vector.h"
#include"showout.h"

int Output_QR(Dmatrix *pQ, Dmatrix *pR)
{
    Dmatrix *pA;
    Dmatrix *pQQT;
    Dmatrix *pQTQ;
    pA = CreateDmatrix(pQ->nDimRow, pR->nDimCol);
    pQTQ = CreateDmatrix(pQ->nDimCol, pQ->nDimCol);
    pQQT = CreateDmatrix(pQ->nDimRow, pQ->nDimRow);

    printf("Matrix Q: \n");
    ShowDmatrix(pQ);
    printf("Matrix R: \n");
    ShowDmatrix(pR);
    
    printf("Matrix QxR: \n");
    DmatrixMulti_NN(pA, pQ, pR);
    ShowDmatrix(pA);

    printf("Matrix Q*xQ: \n");
    DmatrixMulti_TN(pQTQ, pQ, pQ);
    ShowDmatrix(pQTQ);

    printf("Matrix QxQ*: \n");
    DmatrixMulti_NT(pQQT, pQ, pQ);
    ShowDmatrix(pQQT);

    DeleteDmatrix(&pA);
    DeleteDmatrix(&pQQT);
    DeleteDmatrix(&pQTQ);
    return 1;
}
