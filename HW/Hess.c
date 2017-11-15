/* Algorithm 26.1
 * Householder Reduction to Hessenberg Form
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"matrix.h"
#include"vector.h"

void ExtractDmatrixa(Dmatrix *pA, Dmatrix *pB)
{
    
}


void Hess(Dmatrix *pH, Dmatrix *pA)
/* Householder Reduction to Hessenberg Form
 * H = Hess(A)
 * matrix A will not be updated in place.
 */
{
    // B = A
    int m;
    m = pA->nDimCol;
    Dmatrix *pB;
    pB = CreateDmatrix(m ,m);
    CopyDmatrix(pB, pA);
    int k;
    for(k=0; k<m-2; k++)
    {
        // x = B(k+1:m-1, k)
        Dvector *pX, *pV;
        pX = CreateDvector(m-k-1);
        int dim_pX = pX->dim;
        for(int i=0; i < dim_pX; i++)
        {
            pX->data[i] = pB->data[(i+k+1) + k*pB->nDimRow];
        }
        //v_k = sign(x1) ||x|| e1 + x
        double sign;
        sign = pX->data[0] / fabs(pX->data[0]);
        pX->data[0] += sign * Dvector_2Norm(pX);
        pV = pX;
        // v = v/||v||
        DvectorScalar(pV, 1/Dvector_2Norm(pV), pV);
        // update matrix B
    }
}


void main()
{

}
