#!/usr/bin/python
"""
HW7: Algorightm 28.2
"Practical" QR Algorithm
"""

import numpy as np
import sys
import math
from numpy import linalg as LA

def SigExit(**string):
    for i in string:
        print i,
    sys.exit()

def sign(x):
    """return sign(x)"""
    if x >=0:
        return 1
    else:
        return -1

def vvect(v):
    """
    form vetical vector form one-d array
    return 2-d array
    """
    return v.reshape((-1,1))

def hvect(v):
    """
    form horizontal vector from one-d array
    return 2-d array
    """
    return v.reshape((1,-1))

def Hess(matrix_A):
    """
    Household Reduction to Hessenberg From.

    input:  <type: 2d-np.array> matrix_A real or complex
    return: <type: 2d-np.array> H, Q
    Q * H * Q.T = matrix_A
    """
    A = matrix_A.copy().astype(float)
    shape = A.shape
    if shape[0] != shape[1]:
        SigExit("Matrix is not square.")
    dim = shape[0]
    Q = np.identity(dim)
    for k in range(dim-2):
        v_k = A[k+1:, k].copy()
        v_k[0] += sign(v_k[0]) * np.linalg.norm(v_k)
        v_k = v_k/(np.linalg.norm(v_k))
        v_k = vvect(v_k)
        A[k+1:, k:] -= 2*v_k.dot( (v_k.conj().T).dot(A[k+1:, k:]) )
        A[:, k+1:]  -= 2*(A[:, k+1:].dot(v_k)).dot(v_k.conj().T)
        Q[:, k+1:]  -= 2*(Q[:, k+1:].dot(v_k)).dot(v_k.conj().T)
    H = A
    return H, Q

def GQR_Hess(matrix_A):
    """
    QR factorization implemented with givens rotation for matrix with
    Hessenberg form.

    input:  <type: 2d-np.array> matrix_A real or complex
    return: <type: 2d-np.array> Q, R
    Q * R = A
    """
    A = matrix_A.copy()
    shape = A.shape
    if shape[0] != shape[1]:
        SigExit("matrix is no squred.")
    dim = shape[0]
    R = A
    Q = np.identity(dim)
    for i in range(dim-1):
        d = math.sqrt(R[i,i]**2 + R[i+1, i]**2)
        c = R[i,i] / d
        s = R[i+1, i] / d
        J = np.array([[c,s],[-s,c]])
        R[[i, i+1],:] = J.dot(R[[i,i+1],:])
        Q[:,[i,i+1]] = Q[:,[i,i+1]].dot(J.conj().T)
    return Q, R

def Eig_QR(matrix_A):
    """
    "Practical" QR algorithm for eigenvalue decomposition

    input:  <type: 2d-np.array> matrix_A real
    return: <type: 2d-np.array> Sigma, X
    X * Sigma * X.T = matrix_A
    """
    shape = matrix_A.shape
    if shape[0] != shape[1]:
        SigExit("matrix is not square.")
    dim = shape[0]
    A, Q_0 = Hess(matrix_A)
    count = 0
    X = np.identity(dim)
    while True:
        count += 1
        miu = A[dim-1, dim-1]
        for i in range(dim):
            A[i,i] -= miu
        Q, R = GQR_Hess(A)
        A = R.dot(Q)
        for i in range(dim):
            A[i,i] += miu
        X = X.dot(Q)
        if (abs(A[1,0]/A[0,0]) < 10e-15) or (count > 1000):
            print "iter = %d" %count
            X = Q_0.dot(X)
            Sigma = A
            return Sigma, X

def error(A, Sigma, X):
    my_eiw = Sigma.copy()
    my_eiv = X.copy()
    eiw, eiv = LA.eig(A)
    eiw = np.diag(eiw)
    my_eiw -= eiw
    my_eiv -= eiv
    print "error analysis for the whole decomposition:"
    print "error of eigenvalues: %f" %(LA.norm(my_eiw, ord='fro')/LA.norm(eiw, ord='fro'))
    print "error of eigenvectors: %f" %(LA.norm(my_eiv, ord='fro')/LA.norm(eiv, ord='fro'))
    print "error analysis of the first eigenpairs:"
    print "error of eigenvalue: %f" %(my_eiw[0,0]/eiw[0,0])
    print "error of eigenvector: %f" %(LA.norm(my_eiw[:,0])/LA.norm(eiw[:,0]))

def main():
    A = np.random.rand(500,500)
    A = (A+A.T)/2
    np.set_printoptions(precision=3, suppress=True)
    print "From Eig_QR:"
    print "Input matrix A:"
    print A
    Sigma, X = Eig_QR(A)
    print "Eigenvalue matrix:"
    print Sigma
    print "Eigenvector matrix:"
    print X
    error(A, Sigma, X)

if __name__  == "__main__":
    main()
