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

int do_PP_LU(int m, int n, double *inp)
{
    int P[m];
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
    printf("\nmatrix L\n");
    ShowDmatrix(pL);
    printf("\nmatrix U\n");
    ShowDmatrix(pU);
    printf("\nmatrix P\n");
    ShowDmatrix(pP);
    DeleteDmatrixList(4, &pA, &pL, &pU, &pP);
    return 0;
}

int create_HW_matrix(int n, double *out)
{
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            if(i == j)
                out[i+j*n] = 1;
            else if (i>j)
                out[i+j*n] = -1;
            else if (j==n-1)
                out[i+j*n] = 1;
            else
                out[i+j*n] = 0;
        }
    }
    return 0;
}

int main()
{
    printf("**************************\n");
    double inp1[]={2,4,8,6,
                  1,3,7,7,
                  1,3,9,9,
                  0,1,5,8};
    int m=4;
    int n=4;
    do_PP_LU(m, n, inp1);
    printf("**************************\n");

    //8 * 8
    printf("**************************\n");
    double inp2[64];
    create_HW_matrix(8, inp2);
    do_PP_LU(8,8,inp2);
    printf("**************************\n");
    printf("\n");

    // 100 * 100
    printf("**************************\n");
    m=n=100;
    double inp3[10000];
    create_HW_matrix(100, inp3);
    do_PP_LU(100, 100, inp3);
    printf("**************************\n");
    return 0;
}
