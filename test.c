#include<stdio.h>
#include"matrix.h"
#include"vector.h"
#include"QR.h"


int main()
{
    double A[] = {1,2,3, 0,1,0,0,0,1};
    Dmatrix *pA, *pQ, *pR, *pR_inv;
    pA = CreateDmatrix(3,3);
    pQ = CreateDmatrix(3,3);
    pR = CreateDmatrix(3,3);
    pR_inv = CreateDmatrix(3,3);
    InitDmatrix(pA, A);
    printf("A:\n");
    ShowDmatrix(pA);
    MGS_QR(pA, pQ, pR);
    printf("Q:\n");
    ShowDmatrix(pQ);
    printf("R:\n");
    ShowDmatrix(pR);

    Lapack_Dtrtri(pR_inv, pR);
    printf("R-1\n");
    ShowDmatrix(pR_inv);
    TransposeDmatrix(pQ, pQ);
    printf("Q_T:\n");
    ShowDmatrix(pQ);
    DmatrixMulti_NN(pA, pQ, pR);
    ShowDmatrix(pA);
}
