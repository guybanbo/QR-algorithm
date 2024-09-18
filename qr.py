# -*- coding: utf-8 -*-

import numpy.linalg
import numpy as np
from time import time 

def wilkinson_shift(matrix):
    n = matrix.shape[0]
    if n == 1:
        return matrix[0, 0]
    
    # Extract the bottom-right 2x2 submatrix
    submatrix = matrix[-2:, -2:]
    
    # Compute the eigenvalues of the submatrix
    a, b, c, d = submatrix[0, 0], submatrix[0, 1], submatrix[1, 0], submatrix[1, 1]
    trace = a + d
    det = a * d - b * c
    #Quadratic Formula
    mu1 = trace / 2 + np.sqrt((trace / 2)**2 - det)
    mu2 = trace / 2 - np.sqrt((trace / 2)**2 - det)
    
    # Choose the eigenvalue that is closer to the bottom-right element
    if abs(d - mu1) < abs(d - mu2):
        return mu1
    else:
        return mu2
    
    #computes norm of a vector
def vec_norm(v):
    return np.sqrt(np.sum(np.power(v,2)))

def v_vec_norm(v):
    if (v.shape[0] > 1):
        return np.sqrt(np.power(v[0], 2) + np.power(v[1], 2))
    else:
        return np.sqrt(np.power(v[0],2))
    
   
def v_vec_dot(v):
    if (v.shape[0] > 1):
        return np.power(v[0], 2) + np.power(v[1], 2)
    else:
        return v[0]*v[0]
    
    #helps to converge faster, zeros out some elements
def reflection_vector(x):
    v = np.copy(x)
    v[0] += np.sign(x[0]) *v_vec_norm(x)
    eps = 1e-10  # Small epsilon to avoid division by zero
    return v, 2.0 / (v_vec_dot(v) + eps)

def qr_householder(A):
    n = A.shape[0]
    Q = np.eye(n)
    R = np.copy(A)
    
    for j in range(n-1):
        v, c = reflection_vector(R[j:, j])
        
        # Restrict v to the first 2 elements
        v = v[:2]
        
        # Update R
        R_sub = R[j:j+2, j:]
        R[j:j+2, j:] -= c * np.outer(v, v @ R_sub)
        
        # Update Q
        Q_sub = Q[:, j:j+2]
        Q[:, j:j+2] -= c * np.outer(Q_sub @ v, v)

    return Q, R

     
def hessenberg_form(A, epsilon=1e-6):
    
    n = A.shape[0]
    H = np.copy(A)
    Q = np.eye(n)

    for k in range(n - 2):
        x = H[k+1:, k]#extracts the elements in the k-th column under the main diagonal
        v = np.copy(x)
        v[0] += np.sign(x[0]) * vec_norm(x)
        vNorm=vec_norm(v)
        if(vNorm!=0):#added conditon to avoid division by zero
         v /= vec_norm(v)#normalize
        H[k+1:, k:] -= 2.0 * np.outer(v, v @ H[k+1:, k:])#update the lower submatrix starting from (k+1,k), an index under main diag
        H[:, k+1:] -= 2.0 * np.outer(H[:, k+1:] @ v, v)#update columns starting for the k+1 col
        Q[:, k+1:] -= 2.0 * np.outer(Q[:, k+1:] @ v, v)

    # Create a mask for  elements near the diagonal, assign true to the main diagonal and the diagonals above and below it
    mask = np.eye(n, dtype=bool) | np.eye(n, k=1, dtype=bool) | np.eye(n, k=-1, dtype=bool)


    #zeros elements with a false mask value that their absolute values are smaller than epsilon
    H[~mask] = np.where(np.abs(H[~mask]) < epsilon, 0, H[~mask])
    
    return H, Q

def myEigenRecursive(A, epsilon=1e-6):
    n = A.shape[0]
    I = np.eye(n)   
    Q = np.eye(n) 
    
    if n == 1: 
        return A[:, 0], I

    eigenvectors = np.eye(n)

    while np.min(np.abs(np.diag(A, 1))) > epsilon:#continue while the min elemnet abs value in the upper diagonal is bigger than epsilon 
        u = wilkinson_shift(A)# Use Wilkinson's shift
        uI = u * I
        Q, R = qr_householder(A - uI)# Apply QR factorization on the matrix A - uI
        A = R @ Q + uI                  
        eigenvectors = eigenvectors @ Q # update eigenvectors for the final result

    diag_arr_position = np.argmin(np.abs(np.diag(A, 1)))    

    upper_mat = A[:diag_arr_position + 1, :diag_arr_position + 1]
    low_mat = A[diag_arr_position + 1:, diag_arr_position + 1:]

    eigenvalues_upper, eigenvector_upper = myEigenRecursive(upper_mat, epsilon)
    eigenvalues_lower, eigenvector_lower = myEigenRecursive(low_mat, epsilon)
    
    eigenvalues = np.concatenate([eigenvalues_upper, eigenvalues_lower])#concat the eigenvalues returned from the two calls
    
    v1 = np.r_[eigenvector_upper, np.zeros((eigenvector_lower.shape[0], eigenvector_upper.shape[1]))]#vertical concat, padds eigenvector upper with 0's'
    v2 = np.r_[np.zeros((eigenvector_upper.shape[0], eigenvector_lower.shape[1])), eigenvector_lower]
    
    return eigenvalues, eigenvectors @ np.c_[v1, v2]#concat v1 v2 horizontally



def sort_by_same_order(eigenvalues, eigenvectors):
    sorted_indices = np.argsort(eigenvalues)
    sorted_eigenvalues = eigenvalues[sorted_indices]
    sorted_eigenvectors = eigenvectors[:, sorted_indices]
    return sorted_eigenvalues, sorted_eigenvectors


def is_symmetric(matrix):
    matrix = np.array(matrix)
    return np.array_equal(matrix, matrix.T)

if __name__ == '__main__':
    matrix = np.loadtxt("inv_matrix(800 x 800).txt", dtype='f', delimiter=' ')
    print(is_symmetric(matrix))
    print("Converting matrix to Hessenberg form...")
    tt = time()
    hessen_mat, H = hessenberg_form(matrix,1e-6)
    print("Calculating eigenvalues & eigenvectors of matrix using QR Algorithm...") 
    values, vectors = myEigenRecursive(hessen_mat, epsilon=1e-6)
    vectors = H @ vectors
    duration=time()-tt
    print("\nmy_eigen_recursive:\t", duration)
    tt = time()
    eigenvalues, eigenvectors = np.linalg.eig(matrix)
    duration=time()-tt
    print("np.linalg.eig: \t\t", duration)
    sorted_values, sorted_vectors = sort_by_same_order(values, vectors)
    sorted_eigenvalues, sorted_eigenvectors = sort_by_same_order(eigenvalues, eigenvectors)
    np.savetxt('sorted_vectors_python.txt', (np.abs(sorted_vectors)-np.abs(sorted_eigenvectors)), fmt='%.8f', delimiter=' ', comments='') 
    np.savetxt('sorted_values_python.txt', (np.abs(sorted_values)-np.abs(sorted_eigenvalues)), fmt='%.8f', delimiter=' ', comments='')
    print("\nEigenvector norm difference : \t", np.linalg.norm(np.abs(sorted_vectors) - np.abs(sorted_eigenvectors)))
    print("Eigenvalue norm difference  : \t", np.linalg.norm(np.abs(sorted_values) - np.abs(sorted_eigenvalues))) 

