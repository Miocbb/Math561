#!/usr/bin/python
"""
HW7: Algorigthm 26.1
Householder Reduction to Hessenberg Form.
"""
import numpy as np
import sys


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

def main():
    A=np.random.rand(10,10)
    A = (A+A.T)/2
    np.set_printoptions(precision=3, suppress=True)
    print "Input matrix A:"
    print A
    print "\nHessenberg form of matrix A:"
    H, Q = Hess(A)
    print H


if __name__ == "__main__":
    main()
