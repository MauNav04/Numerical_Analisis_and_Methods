import numpy as np

def angulo(A, i, j):
    Aij = A[i, j]
    Aii = A[i, i]
    Ajj = A[j, j]
    if abs(Aij) > 1e-16:
        theta = 0.5 * np.arctan(2 * Aij / (Aii - Ajj))
    else:
        theta = 0
    return theta

def matrix_rotation(i, j, n, theta):
    G = np.eye(n)
    G[i, i] = np.cos(theta)
    G[j, j] = np.cos(theta)
    G[i, j] = -np.sin(theta)
    G[j, i] = np.sin(theta)
    return G

def jacobi_valores_propios(A, iterMax, tol):
    n = A.shape[0]
    Ak = A.copy().astype(float)
    xk = np.diag(Ak).copy()   # valores propios del paso anterior
    ek = np.inf

    for k in range(iterMax):
        # Encontrar el mayor elemento fuera de la diagonal
        Ak_off = np.abs(Ak - np.diag(np.diagonal(Ak)))
        i, j = np.unravel_index(np.argmax(Ak_off), Ak.shape)

        # Criterio de parada basado en el elemento fuera de diagonal
        if abs(Ak[i, j]) < tol:
            break

        theta = angulo(Ak, i, j)
        G = matrix_rotation(i, j, n, theta)
        Ak = G.T @ Ak @ G

        xk_new = np.diag(Ak).copy()

        # ✅ Error correcto: norma 1 de la diferencia del vector de valores propios
        ek = np.linalg.norm(xk_new - xk, ord=1)

        xk = xk_new

        if ek < tol:
            break

    return xk, ek

# Matriz A 15x15 con A_ij = 0.5*(i+j), i,j = 1,...,15
n = 15
A = np.fromfunction(lambda i, j: 0.5 * (i + j + 2), (n, n))

iterMax = 10000
tol = 1e-10

xk, ek = jacobi_valores_propios(A, iterMax, tol)

print("Valores propios aproximados:")
print(np.sort(xk)[::-1])
print(f"\nError final: {ek:.16e}")