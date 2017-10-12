#include<stdio.h>
#include<stdlib.h>

#include"vector.h"
#include"matrix.h"
#include"showout.h"

int MGS_QR(Dmatrix *pA, Dmatrix *pQ, Dmatrix *pR)
{
    /* Start Modified Gram-Schmidt QR factorization*/
    Dmatrix *pV;
    pV = CreateDmatrix(pA->nDimRow, pA->nDimCol);
    Dvector *pA_i, *pQ_i, *pV_i;
    pA_i = CreateDvector(pA->nDimRow);
    pQ_i = CreateDvector(pQ->nDimRow);
    pV_i = CreateDvector(pA->nDimRow);

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
        for(j=i+1; j < pA->nDimCol; j++)
        {
            ExtractDmatrixCol(pV_i, pV, j);
            r_ij = DvectorIP(pQ_i, pV_i);
            pR->data[i + j*(pR->nDimRow)] = r_ij;
            DvectorArithmetic(pV_i, 1, pV_i, -r_ij, pQ_i);
            InitDmatrixCol(pV, pV_i->data, j);
        }
    }
     /*DO full QR factorization*/
    /*
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
    */
    DeleteDvector(&pA_i);
    DeleteDvector(&pQ_i);
    DeleteDvector(&pV_i);
    return 1;
}
