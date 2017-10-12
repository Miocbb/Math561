#include"least_square.h"
#include"matrix.h"
#include"vector.h"
#include"QR.h"

void LeastSquare(Dmatrix *pA, Dvector *pX, Dvector *pY)
{
    Dmatrix *pQ, *pR, *pQ_T; 
    Dvector *pQY;
    
    pQ = CreateDmatrix(pA->nDimRow, pA->nDimCol);
    pR = CreateDmatrix(pA->nDimCol, pA->nDimCol);
    pQ_T = CreateDmatrix(pA->nDimCol, pA->nDimRow);
    pQY = CreateDvector(pA->nDimCol);

    MGS_QR(pA, pQ, pR);
    TransposeDmatrix(pQ_T, pQ);
    Lapack_Dgemv(pQY, pQ_T, pY); /* QY = Q * Y*/
    SolveUpTriMatrixEq(pR, pX, pQY); /* solve RX = Q_T * Y */

    DeleteDmatrixList(3, &pQ, &pR, &pQ_T);
    DeleteDvectorList(1, &pQY);
}
