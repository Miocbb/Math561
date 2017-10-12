#include<stdio.h>
#include<stdlib.h>

#include"debug.h"
#include"vector.h"
#include"matrix.h"
#include"showout.h"

#define ARRAY_SIZE(array) (int)( sizeof(array)/sizeof(array[0]) )

//#define DEBUG
int GS_QR(double *pInput, int InputRow, int InputCol);

int main()
{
    int InputRow, InputCol;
    double Input1[3][2]={{1,2},
                         {0,1},
                         {1,0}};
    double Input2[3][2]={{1,2},
                         {2,1},
                         {2,2}};
    
    printf("First Input\n");
    InputRow = ARRAY_SIZE(Input1);
    InputCol = ARRAY_SIZE(Input1[0]);
    GS_QR(Input1[0], InputRow, InputCol);
    
    printf("\n\n");
    printf("Second Input\n");
    InputRow = ARRAY_SIZE(Input2);
    InputCol = ARRAY_SIZE(Input2[0]);
    GS_QR(Input2[0], InputRow, InputCol);
    return 1;
}


int GS_QR(double *pInput, int InputRow, int InputCol)
{
    #ifdef DEBUG
        debug();
    #endif
    
    printf("*******************************\n");
    printf("   Classical Gram-Schmidt QR   \n");
    printf("*******************************\n");
    
    Dmatrix *pA, *pQ, *pR;

    pA = CreateDmatrix(InputRow, InputCol);
    pQ = CreateDmatrix(InputRow, InputCol);
    pR = CreateDmatrix(InputCol, InputCol);
    
    //InitDmatrix(pA, &a[0][0]);
    InitDmatrix(pA,pInput);
    printf("Input Matrix A: \n");
    ShowDmatrix(pA);
    
    /* Start classic Gram-Schmidt QR factorization*/
    Dvector *pA_i, *pQ_i, *pV_i;
    pA_i = CreateDvector(pA->nDimRow);
    pQ_i = CreateDvector(pQ->nDimRow);
    pV_i = CreateDvector(pA->nDimRow);

    /* reduced GS QR factorization */
    int i, j;
    double r_jj; /* the diagonal element of R*/
    double r_ij;
    for(j=0; j<InputCol; j++)
    {
        ExtractDmatrixCol(pA_i, pA, j);
        CopyDvector(pV_i, pA_i);
        
        if(j > 0){
            for(i=0; i<j; i++)
            {
                ExtractDmatrixCol(pQ_i, pQ, i);
                r_ij = DvectorIP(pQ_i, pA_i);
                pR->data[j*(pR->nDimRow)+i] = r_ij;
                DvectorArithmetic(pV_i, 1, pV_i, -r_ij, pQ_i);
            }
        }
        r_jj = Dvector_2Norm(pV_i);
        pR->data[j*(pR->nDimCol)+j] = r_jj;

        /* taking low-rank input matrix into consideration*/
        if(r_jj != 0){
            DvectorScalar(pV_i, 1/r_jj, pV_i);
        }else{
            int k;
            double norm = 0;
            RandomDvector(pV_i);
            while(norm == 0)
            {
                for(k=0; k<j; k++)
                {   
                    ExtractDmatrixCol(pQ_i, pQ, k);
                    DvectorArithmetic(pV_i, 1, pV_i, -DvectorIP(pQ_i, pV_i), pQ_i);
                }
                norm = Dvector_2Norm(pV_i);
            }
            DvectorScalar(pV_i, 1/norm, pV_i);
        }
        InitDmatrixCol(pQ, pV_i->data, j);
    }


    /*Show reduced QR factorization*/
    printf("****************************** \n");
    printf("Show reduced QR factorization \n");
    printf("****************************** \n");
    Output_QR(pQ, pR);

    /*DO full QR factorization*/
    printf("****************************** \n");
    printf("  Show full QR factorization\n");
    printf("****************************** \n");
    if(pA->nDimCol<pA->nDimRow)
    {
        Dmatrix *pQ_Expd; 
        Dvector *pQ_Expd_i;
        double norm;
        int i, col_Q_Expd;
        col_Q_Expd = pQ->nDimRow - pQ->nDimCol;

        pQ_Expd = CreateDmatrix(pQ->nDimRow, col_Q_Expd);
        pQ_Expd_i  = CreateDvector(pQ->nDimRow);
        for(i=0; i< col_Q_Expd; i++)
        {
            int k, count=0;
            norm = 0;
            while(norm < 0.0001)
            {
                count++;
                RandomDvector(pQ_Expd_i);
                for(k=0; k < pQ->nDimCol; k++)
                {
                    ExtractDmatrixCol(pQ_i, pQ, k);
                    DvectorArithmetic(pQ_Expd_i, 1, pQ_Expd_i, -DvectorIP(pQ_i, pQ_Expd_i), pQ_i);
                }
                norm = Dvector_2Norm(pQ_Expd_i);
            }
            printf("Cycle %d times to find othogonal vector to do full QR!\n", count);
            DvectorScalar(pQ_Expd_i, 1/norm, pQ_Expd_i);
            InitDmatrixCol(pQ_Expd, &(pQ_Expd_i->data[0]), i);
        }
        DmatrixExpansionByCol(pQ, pQ_Expd);
        DeleteDmatrix(&pQ_Expd);

        Dmatrix *pR_Expd;
        pR_Expd = CreateDmatrix(pA->nDimRow-pA->nDimCol, pA->nDimCol);
        DmatrixExpansionByRow(pR, pR_Expd);
        DeleteDmatrix(&pR_Expd);
    }

    DeleteDvector(&pA_i);
    DeleteDvector(&pQ_i);
    DeleteDvector(&pV_i);

    /*Show Output_QR*/
    Output_QR(pQ, pR);
    return 1;
}

