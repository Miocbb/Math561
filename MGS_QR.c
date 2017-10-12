#include<stdio.h>
#include<stdlib.h>

#include"vector.h"
#include"matrix.h"
#include"showout.h"

#define ARRAY_SIZE(array) (int)( sizeof(array)/sizeof(array[0]) )

//#define DEBUG
int MGS_qr(double *pInput, int InputRow, int InputCol);

int main()
{
    double Input1[]={1,0,1,2,1,0};
    double Input2[]={1,2,2,1,2,2};
    
    printf("First Input\n");
    MGS_qr(Input1, 3,2);
    
    printf("\n\n");
    printf("Second Input\n");
    MGS_qr(Input2,3,2);
    return 1;
}


int MGS_qr(double *pInput, int InputRow, int InputCol)
{
    printf("*******************************\n");
    printf(" Modified GS QR factorization\n");
    printf("*******************************\n");
    
    Dmatrix *pA, *pQ, *pR, *pV;
    pA = CreateDmatrix(InputRow, InputCol);
    pV = CreateDmatrix(InputRow, InputCol);
    pQ = CreateDmatrix(InputRow, InputCol);
    pR = CreateDmatrix(InputCol, InputCol);
    
    InitDmatrix(pA, pInput);
    printf("Input Matrix A: \n");
    ShowDmatrix(pA);
    /* Start Modified Gram-Schmidt QR factorization*/
    Dvector *pA_i, *pQ_i, *pV_i;
    pA_i = CreateDvector(pA->nDimRow);
    pQ_i = CreateDvector(pQ->nDimRow);
    pV_i = CreateDvector(pA->nDimRow);

    printf("****************************** \n");
    printf("Show reduced MGS QR factorization \n");
    printf("****************************** \n");
    /* reduced MGS QR factorization */
    int i, j;
    double r_jj; /* the diagonal element of R*/
    double r_ij;

    CopyDmatrix(pV, pA);
    for(i=0; i < pA->nDimCol; i++)
    {
        ExtractDmatrixCol(pV_i, pV, i);
        r_jj = Dvector_2Norm(pV_i);
        pR->data[i + i*(pR->nDimRow)] = r_jj;
        if(r_jj == 0)
        {
            /* taking low-rank A into consideration
             * generate orthorgonal pQ_i;*/ 
            double norm = 0;
            while(norm <0.00001){
                RandomDvector(pV_i);
                if(i != 0)
                {
                    int k;
                    for(k=0; k<i; k++)
                    {
                        ExtractDmatrixCol(pQ_i, pQ, k);
                        DvectorArithmetic(pV_i, 1, pV_i, -DvectorIP(pQ_i, pV_i), pQ_i);
                    }
                }
                norm = Dvector_2Norm(pV_i);
            }
            DvectorScalar(pQ_i, 1/norm, pV_i);
        }
        else{
            DvectorScalar(pQ_i, 1/r_jj, pV_i);
        }
        InitDmatrixCol(pQ, pQ_i->data, i);
        /*update V column vector: V_i+1 to V_n*/
        for(j=i+1; j < pV->nDimCol; j++)
        {
            ExtractDmatrixCol(pV_i, pV, j);
            r_ij = DvectorIP(pQ_i, pV_i);
            pR->data[i + j*(pR->nDimRow)] = r_ij;
            DvectorArithmetic(pV_i, 1, pV_i, -r_ij, pQ_i);
            InitDmatrixCol(pV, pV_i->data, j);
        }
    }
    Output_QR(pQ, pR);

    /*DO full QR factorization*/
    printf("****************************** \n");
    printf("Show full MGS QR factorization\n");
    printf("****************************** \n");
    if(pA->nDimCol < pA->nDimRow)
    {
        Dmatrix *pQ_Expd; 
        Dvector *pQ_Expd_i;
        double norm;
        int i, col_Q_Expd;
        col_Q_Expd = pQ->nDimRow - pQ->nDimCol;

        pQ_Expd = CreateDmatrix(pQ->nDimRow, 1);
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
            DmatrixExpansionByCol(pQ, pQ_Expd);
        }
        DeleteDmatrix(&pQ_Expd);

        Dmatrix *pR_Expd;
        pR_Expd = CreateDmatrix(pA->nDimRow-pA->nDimCol, pA->nDimCol);
        DmatrixExpansionByRow(pR, pR_Expd);
        DeleteDmatrix(&pR_Expd);
    }
    /*Show Output_QR*/
    Output_QR(pQ, pR);
    DeleteDvector(&pA_i);
    DeleteDvector(&pQ_i);
    DeleteDvector(&pV_i);
    return 1;
}
