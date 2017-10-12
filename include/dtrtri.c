#include<stdio.h>
#include"matrix.h"
#include"vector.h"
#include"lapack.h"

void Lapack_Dtrtri(Dmatrix * pA, Dmatrix *pB)
    /* A = B-1
     * B is an upper triangular matrix
     * */
{   
    if(pA == NULL || pB == NULL)
    {
        printf("Lapack_Dtritir: matrix is NULL, inverse failed!\n");
        return;
    }
    Dmatrix *pB_tem;
    int INFO;
    pB_tem = CreateDmatrix(pB->nDimRow, pB->nDimCol);
    CopyDmatrix(pB_tem, pB);
    dtrtri_("U", "N", &(pB_tem->nDimCol), pB_tem->data, 
            &(pB_tem->nDimCol), &INFO);
    CopyDmatrix(pA, pB_tem);
    DeleteDmatrix(&(pB_tem));
}
