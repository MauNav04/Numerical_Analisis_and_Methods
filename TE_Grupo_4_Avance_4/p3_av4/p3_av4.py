import numpy as np

def metodo_qr(A, tol, max_iter):
    n = A.shape[0]
    Ak = A.copy().astype(float)
    xk = np.diag(Ak).copy()   # vector de valores propios anterior
    ek = np.inf
    iters = 0

    for k in range(max_iter):
        Q, R = np.linalg.qr(Ak)
        Ak = R @ Q

        xk_new = np.diag(Ak).copy()

        # ✅ Error correcto: norma 1 de la diferencia del VECTOR de valores propios
        ek = np.linalg.norm(xk_new - xk, ord=1)

        xk = xk_new
        iters += 1

        if ek < tol:
            break

    return xk, ek, iters

# Matriz A 20x20
n = 20
A = np.zeros((n, n))
for i in range(n):
    for j in range(n):
        if i == j:
            A[i, j] = 2
        elif abs(i - j) == 1:
            A[i, j] = -1

tol      = 1e-10
max_iter = 10000

valores_propios, error_final, iteraciones = metodo_qr(A, tol, max_iter)

print("Valores propios aproximados:")
print(valores_propios)
print(f"Error final: {error_final:.16e}")
print(f"Número de iteraciones realizadas: {iteraciones}")