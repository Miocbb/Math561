#include<stdio.h>

#include"vector.h"
#include"matrix.h"
#include"debug.h"

#define ARRAY_SIZE(array) ( sizeof(array)/sizeof(array[0]) )

void debug()
{
    double a[]={1,2,3};
    double b[]={2,1,2};

    Dvector *pa, *pb, *pc;
    
    printf("size of a and b are %d, %d \n",(int)sizeof(a), (int)sizeof(b));
    pa = CreateDvector(ARRAY_SIZE(a));
    pb = CreateDvector(ARRAY_SIZE(b));
    pc = CreateDvector(ARRAY_SIZE(b));
    
    InitDvector(pa, a);
    InitDvector(pb, b);
    printf("Below is vector a: \n");
    ShowDvector(pa);
    printf("Below is vector b: \n");
    ShowDvector(pb);
    
    printf("Below is vector a+b: \n");
    DvectorArithmetic(pc, 1, pa, 1, pb);
    ShowDvector(pc);
    printf("Below is vector a-b: \n");
    DvectorArithmetic(pc, 1, pa, -1, pb);
    ShowDvector(pc);
    printf("Below is vector a*2: \n");
    DvectorScalar(pc, 2, pa);
    ShowDvector(pc);
    printf("Below is vector a and b 2-norm: \n");
    printf("%5f,  %5f\n", Dvector_2Norm(pa), Dvector_2Norm(pb));

    double A[]={1,2,3,1,0,2};
    Dmatrix *pA;
    Dvector *pX;
    
    pA = CreateDmatrix(2,3);
    pX = CreateDvector(pA->nDimRow);
    InitDmatrix(pA, A);
    printf("size of matrix A: Row=%d, Col=%d\n", pA->nDimRow, pA->nDimCol);
    
    printf("Below is Matrix A:\n");
    ShowDmatrix(pA);
    
    printf("Below is 1-col vector of A:\n");
    ExtractDmatrixCol(pX, pA, 0);
    ShowDvector(pX);
    DeleteDvector(&pX);

    printf("Below is 1-row vector of A:\n");
    pX = CreateDvector(pA->nDimCol);
    ExtractDmatrixRow(pX, pA, 0);
    ShowDvector(pX);
    DeleteDvector(&pX);

    double rowVect[]={111, 222,333};
    double colVect[]={444,555};
    printf("initial matrix 1-row\n");
    InitDmatrixRow(pA, rowVect, 0);
    ShowDmatrix(pA);
    
    printf("initial matrix 1-col\n");
    InitDmatrixCol(pA, colVect, 0);
    ShowDmatrix(pA);

    double B[]={1,1,2,2};
    double C[]={1,2,1,2,1,2};
    Dmatrix *pB, *pC;
    pB = CreateDmatrix(2,2);
    pC = CreateDmatrix(2,3);
    InitDmatrix(pA, A);
    InitDmatrix(pB, B);
    InitDmatrix(pC, C);
    DmatrixExpansionByCol(pA, pB);
    ShowDmatrix(pA);
    DeleteDmatrix(&pA);
    pA = CreateDmatrix(2,3);
    InitDmatrix(pA, A);
    printf("A+C expansion\n");
    ShowDmatrix(pA);
    ShowDmatrix(pC);
    DmatrixExpansionByRow(pA, pC);
    printf("\n");
    ShowDmatrix(pA);
    
    Dmatrix *pAT;
    pAT=CreateDmatrix(pA->nDimCol, pA->nDimRow);
    TransposeDmatrix(pAT, pA);
    ShowDmatrix(pAT);

    double A1[]={1,0,1,1,2,3};
    double B1[]={1,2,0,3,0,0,0,0,0};
    DeleteDmatrix(&pA);
    DeleteDmatrix(&pB);
    DeleteDmatrix(&pC);
    pA=CreateDmatrix(2,3);
    pB=CreateDmatrix(3,3);
    pC=CreateDmatrix(3,3);
    InitDmatrix(pA, A1);
    InitDmatrix(pB, B1);
    printf("matrix A1\n");
    ShowDmatrix(pA);
    printf("matrix B1\n");
    ShowDmatrix(pB);
    
    DmatrixMulti_TN(pC, pB, pB);
    printf("A1*B1\n");
    ShowDmatrix(pC);

    DeleteDmatrix(&pC);
    ShowDmatrix(pC);
}
