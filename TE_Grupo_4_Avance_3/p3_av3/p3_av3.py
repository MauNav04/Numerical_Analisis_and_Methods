import numpy as np

def factLU(A):
    A = np.array(A, dtype=float)
    n = A.shape[0]
    L = np.eye(n)
    U = A.copy() 
    # Realizar la factorización LU por Gauss
    for k in range(n - 1):
        for i in range(k + 1, n):
            L[i, k] = U[i, k] / U[k, k]
            U[i, k:] -= L[i, k] * U[k, k:]
    return L, U

def sustitucion_adelante(L, b):
    n = len(b)
    y = np.zeros(n)
    for i in range(n):
        y[i] = b[i] - np.dot(L[i, :i], y[:i]) # L*y=b
    return y

def sustitucion_atras(U, y):
    n = len(y)
    x = np.zeros(n)
    # U*x=y
    for i in range(n - 1, -1, -1):
        x[i] = (y[i] - np.dot(U[i, i+1:], x[i+1:])) / U[i, i]
    return x

def resolver_LU(L, U, b):
    return sustitucion_atras(U, sustitucion_adelante(L, b))

# Construir matriz tridiagonal 1000x1000
n = 1000
A = np.diag(5*np.ones(n)) + np.diag(np.ones(n-1), 1) + np.diag(np.ones(n-1), -1)

# Factorización LU (una sola vez)
L, U = factLU(A)

# Simulación: A x^(k+1) = x^(k)  =>  x^(k+1) = A^(-1) x^(k)
x = np.ones(n)
for k in range(1, 10001):
    x = resolver_LU(L, U, x)
    norma = np.linalg.norm(x)
    print(f"k={k}, ||x||={norma:.4e}")
    # Se detiene si norma<tol
    if norma < 1e-8:
        print(f"\nConvergencia en iteración {k}: ||x||={norma:.4e} < 1e-8")
        print("Sistema estable y disipativo verificado.")
        break