/* This is the implementation of LU factorization with partial
 * pivoting.
 */

#include<stdio.h>
#include<stdlib.h>
#include"matrix.h"
#include"vector.h"

void MakePivotingMatrix(int *P, Dmatrix *pP)
    /* convert array stored pivoting matrix P
     * back to normal dense matrix, so that implemented
     * function "ShowDmatrix" can be directly used to
     * show output.
     * */
{
    int i, j;
    int dim=pP->nDimCol;
    for(i=0; i<dim; i++)
    {
        for(j=0; j<dim; j++)
        {
            if(j+1 == P[i])
                pP->data[j*dim+i] = 1;
            else
                pP->data[j*dim+i] = 0;
        }
    }
}


void MakeIdentityDMatrix(Dmatrix *pA)
    /* A =I */
{
    int dim=pA->nDimCol;
    int i,j;
    for(i=0; i<dim; i++)
        for(j=0; j<dim; j++)
        {
            pA->data[i*dim+j]=0;
            pA->data[i*dim+i]=1;
        }
}

void PP_LU(Dmatrix *pA, Dmatrix *pL, Dmatrix *pU, int *pP)
    /* LU factorization with partial pivoting:
     * PA = LU
     * note: pivoting matrix is stored by int array with
     * row elementary unit vectors {e_i}.
     * e.g. P = [1,2,3,4] = [e1, e2, e3, e4]^T
     */
{
    int k;
    int m, n;
    n = pA->nDimCol;
    m = pA->nDimRow;
    /* initialization*/
    CopyDmatrix(pU, pA); //U=A
    MakeIdentityDMatrix(pL); //L=I
    for(int count=0; count<m; count++)//P=I
        pP[count] = count+1;

    for(k=0; k < m-1; k++)
    {
        // i = argmax(Ujk), where k<=j<=m-1
        int i, j;
        i = k;
        for(j=k; j<m; j++)
        {
            if(pU->data[k*m+i] < pU->data[k*m+j])
                i = j;
        }

        // U(k, k:n-1) <=> U(i, k:n-1)
        double tem_double;
        for(int count=k; count<n; count++)
        {
            tem_double = pU->data[count*m+k];
            pU->data[count*m+k] = pU->data[count*m+i];
            pU->data[count*m+i] = tem_double;
        }

        // if k>=0:
        // L(k, 0:k-1) <=> L(i, 0:k-1)
        if(k>0)
        {
            for(int count=0; count<k; count++)
            {
                tem_double = pL->data[count*m+k];
                pL->data[count*m+k] = pL->data[count*m+i];
                pL->data[count*m+i] = tem_double;
            }
        }

        // P[k] <=> P[i]
        int tem_int;
        tem_int = pP[k];
        pP[k] = pP[i];
        pP[i] = tem_int;

        // Update other rows
        for(int j=k+1; j<m; j++)
        {
            //L(j,k) = U(j,k)/U(k,k)
            pL->data[k*m+j] = pU->data[k*m+j]/pU->data[k*m+k];
            //U(j, k:n-1) = U(j,k:n-1) - L(j,k)U(k:n-1)
            for(int count=k; count<n; count++)
            {
                pU->data[count*m+j] -= pL->data[k*m+j] * pU->data[count*m+k];
            }
        }
    }
}


int main()
{
    double inp[]={1,2,3,
                  4,5,6,
                  7,8,9,
                  10,11,12};
    int m=3;
    int n=4;
    int P[3];
    Dmatrix *pA, *pL, *pU, *pP;

    pA = CreateDmatrix(m,n);
    pL = CreateDmatrix(m,m);
    pU = CreateDmatrix(m,n);
    pP = CreateDmatrix(m,m);
    InitDmatrix(pA, inp);

    PP_LU(pA, pL, pU, P);
    MakePivotingMatrix(P, pP);
    printf("PA = LU\n");
    printf("matrix A\n");
    ShowDmatrix(pA);
    printf("matrix L\n");
    ShowDmatrix(pL);
    printf("matrix U\n");
    ShowDmatrix(pU);
    printf("matrix P\n");
    ShowDmatrix(pP);
    DeleteDmatrixList(4, &pA, &pL, &pU, &pP);
    return 0;
}
