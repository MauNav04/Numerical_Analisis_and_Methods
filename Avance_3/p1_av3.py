import numpy as np
import math
import time
import matplotlib.pyplot as plt

# Función para construir la matriz A
def construir_A(m=45, n=30):
    i = np.arange(1, m + 1).reshape(-1, 1)
    j = np.arange(1, n + 1).reshape(1, -1)
    return i**2 + j**2

# Función para calcular la pseudoinversa de Li y Li
def pseudoinversa_li_li(A, p, tol=1e-5, maxit=10000):
    m, n = A.shape
    # Valor inicial tipo Newton-Schulz
    alpha = 1.0 / (np.linalg.norm(A, 1) * np.linalg.norm(A, np.inf))
    X = alpha * A.T

    err = np.linalg.norm(A - A @ X @ A, ord='fro')
    k = 0

    t0 = time.perf_counter()

    while err > tol and k < maxit:
        M = A @ X  # m x m

        S = np.zeros((m, m))
        I = np.eye(m)

        # Para q=1, la potencia es M^0 = I
        M_power = I.copy()

        for q in range(1, p + 1):
            coef = ((-1) ** (q - 1)) * math.comb(p, q)
            if q == 1:
                M_power = I
            elif q == 2:
                M_power = M
            else:
                M_power = M_power @ M

            S += coef * M_power

        X_new = X @ S

        err = np.linalg.norm(A - A @ X_new @ A, ord='fro')

        X = X_new
        k += 1

    tf = time.perf_counter()

    return X, err, k, tf - t0

# Función para resolver mínimos cuadrados
def resolver_minimos_cuadrados(X, A, b):
    x = X @ b
    residual = np.linalg.norm(A @ x - b)
    return x, residual

# Función para graficar los resultados
def graficar_resultados(resultados):
    p_vals = [r["p"] for r in resultados]
    it_vals = [r["iteraciones"] for r in resultados]
    t_vals = [r["tiempo"] for r in resultados]

    plt.figure()
    plt.plot(p_vals, it_vals, marker='o')
    plt.xlabel("p")
    plt.ylabel("Número de iteraciones")
    plt.title("p vs número de iteraciones")
    plt.grid(True)
    plt.savefig("p_vs_iteraciones.png", dpi=300, bbox_inches='tight')

    plt.figure()
    plt.plot(p_vals, t_vals, marker='o')
    plt.xlabel("p")
    plt.ylabel("Tiempo de ejecución (s)")
    plt.title("p vs tiempo de ejecución")
    plt.grid(True)
    plt.savefig("p_vs_tiempo.png", dpi=300, bbox_inches='tight')

    plt.show()

# Función para imprimir los resultados en una tabla
def imprimir_tabla(resultados):
    print("p\tError Final\tIteraciones\tTiempo(s)\tResidual")
    for r in resultados:
        print(f"{r['p']}\t{r['error_final']:.6e}\t{r['iteraciones']}\t{r['tiempo']:.6f}\t{r['residual']:.6e}")

# Función principal que corre el experimento
def experimento():
    A = construir_A()
    b = np.ones((45, 1))

    valores_p = [1, 2, 3, 4, 5, 6, 7, 8, 10]
    resultados = []

    for p in valores_p:
        X, err_final, iters, tiempo = pseudoinversa_li_li(A, p, tol=1e-5, maxit=10000)
        x, residual = resolver_minimos_cuadrados(X, A, b)

        resultados.append({
            "p": p,
            "X": X,
            "error_final": err_final,
            "iteraciones": iters,
            "tiempo": tiempo,
            "x": x,
            "residual": residual
        })

        print(f"p={p} | error={err_final:.6e} | it={iters} | tiempo={tiempo:.6f} s | residual={residual:.6e}")

    imprimir_tabla(resultados)
    graficar_resultados(resultados)

if __name__ == "__main__":
    experimento()