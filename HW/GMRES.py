#!/usr/bin/python
"""
HW8: Algorithm 35.1
Implementation of GMRES Algorithm.
"""

import sys
import math
import numpy as np

def SigExit(*string):
    for i in string:
        print i,
    sys.exit()

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

def GQR_Hess(mat_A):
    """
    QR factorization implemented with givens rotation.
    -matrix is upper Hessenberg form.
    -reduced QR is returned.

    input:  <type: 2d-np.array> matrix_A real or complex
    return: <type: 2d-np.array> Q, R
    Q * R = A
    """
    A = mat_A.copy().astype(float)
    m = A.shape[0]
    n = A.shape[1]
    R = A
    Q = np.identity(m)
    for i in range(n):
        if i+1 <= m-1:
            d = math.sqrt(R[i,i]**2 + R[i+1, i]**2)
            c = R[i,i] / d
            s = R[i+1, i] / d
            J = np.array([[c,s],[-s,c]])
            R[[i, i+1], i:] = J.dot(R[[i,i+1], i:])
            Q[:,[i,i+1]] = Q[:,[i,i+1]].dot(J.conj().T)
    return Q[:, :n], R[:n, :]


def GMRES(mat_A, vec_b):
    """
    solve linear equation with GMRES method.
    A is square matrix.

    input: <type: 2d-np.array> matrix_A
           <type: 1d-np.arrray> vector_b
    return: <type: 1d-np.array> vector_X
    A * X = b
    """
    A = mat_A.astype(float)
    b = vvect(vec_b).astype(float)
    dim_A = A.shape
    if dim_A[0] != dim_A[1]:
        SigExit("Terminated: not square matrix input.")
    q = list() # store Arnoldi basis for Krylov space
    b_norm = np.linalg.norm(b)
    q.append(b/b_norm)
    Q_n = q[0] # align the Arnoldi basis to be a matrix
    for n in range(dim_A[0]):
        niter = n + 1
        V = A.dot(q[n])
        H = np.zeros((n+2,n+1))
        if n == 0:
            H_old = H.copy()
        if n > 0 :
            H[0:n+1, 0:n] = H_old # update H matrix
        # Arnoldi iteration
        for j in range(n+1):
            H[j, n] = q[j].conj().T.dot(V)
            V -= H[j, n] * q[j]
        H[n+1, n] = np.linalg.norm(V)
        H_old = H.copy()
        if H[n+1, n] <= 1e-14:
            break # Krylov space complete to exit.
        q.append(V / H[n+1, n])
        # start getting approximate solution for AX = b
        Q,R = GQR_Hess(H)
        Y = b_norm * (np.linalg.inv(R).dot(Q.conj().T)[:,0])
        Y = vvect(Y)
        X = Q_n.dot(Y)
        # check residual
        rsd = H.dot(Y)
        rsd[0] -= b_norm
        rsd_norm = np.linalg.norm(rsd)
        if rsd_norm / b_norm  <  1e-10:
            return X[:, 0], niter
        Q_n = np.hstack((Q_n,q[-1])) # update basis matrix of Krylov space
    return X[:,0], niter


def main():
    A = np.zeros((805,805))
    B = np.random.rand(805,805)
    Q,R = np.linalg.qr(B)
    for i in range(801):
        A[i,i] = 1 + i/100
    for i in range(4):
        A[801+i,801+i] = 10 + 2**i
    A = Q.dot(A).dot(Q.T)
    b = np.random.rand(805)
    np.set_printoptions(precision=3, suppress=True)
    X, niter = GMRES(A,b)
    print "niter = ", niter
    print "first 10 elements of solution"
    print X[:10]
    print "first 10 elements of b"
    print b[:10]
    print "first 10 elementf of AX"
    print A.dot(vvect(X))[:,0][:10]

if __name__ == "__main__":
    main()
