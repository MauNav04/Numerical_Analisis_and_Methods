import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Ecuacion de onda unidimensional
# ============================================================

def wave1d(f, g, h, k, T):
    # --------------------------------------------------------
    # PASO 1: Calcular nodos espaciales y temporales
    # --------------------------------------------------------
    m = int(round(1.0 / h))   # numero de subintervalos espaciales
    n = int(round(T / k))     # numero de subintervalos temporales

    x = np.zeros(m + 1)
    for i in range(m + 1):
        x[i] = i * h

    t = np.zeros(n + 1)
    for j in range(n + 1):
        t[j] = j * k

    # --------------------------------------------------------
    # PASO 2: Parametro r = k/h (numero de Courant)
    # --------------------------------------------------------
    r = k / h

    print(f"Parametro de Courant r = k/h = {r:.6f}")
    if r > 1:
        print("ADVERTENCIA: r > 1, el metodo puede ser inestable.")
    else:
        print("Condicion de estabilidad r <= 1: SATISFECHA")

    # --------------------------------------------------------
    # PASO 3: Inicializar matriz U de tamaño (m+1) x (n+1)
    # --------------------------------------------------------
    U = np.zeros((m + 1, n + 1))

    # --------------------------------------------------------
    # PASO 4: Condicion inicial j=0  U[i,0] = f(x_i)
    # --------------------------------------------------------
    for i in range(m + 1):
        U[i, 0] = f(x[i])

    # --------------------------------------------------------
    # PASO 5: Primer paso temporal j=0 → j=1
    # --------------------------------------------------------
    for i in range(1, m):
        U[i, 1] = (U[i, 0]
                   + k * g(x[i])
                   + (r**2 / 2) * (U[i+1, 0] - 2*U[i, 0] + U[i-1, 0]))

    # --------------------------------------------------------
    # PASO 6: Avance en el tiempo j=1, 2, ..., n-1
    # --------------------------------------------------------
    for j in range(1, n):
        for i in range(1, m):
            U[i, j+1] = (2*U[i, j] - U[i, j-1]
                         + r**2 * (U[i+1, j] - 2*U[i, j] + U[i-1, j]))

    return x, t, U


# ============================================================
# APLICACION: resolver el problema dado
# ============================================================

# Condiciones iniciales del problema
def f(x):
    return np.sin(np.pi * x)       # u(x, 0) = sin(pi*x)

def g(x):
    return 0.0                     # du/dt(x, 0) = 0

# Solucion exacta
def u_exacta(x, t):
    return np.cos(np.pi * t) * np.sin(np.pi * x)

# --------------------------------------------------------
# Seleccion de h y k garantizando estabilidad (r = k/h <= 1)
# Se elige h = 0.1  y  k = 0.1  →  r = 1.0  (estable)
# --------------------------------------------------------
h = 0.1
k = 0.1
T = 1.0

print("=" * 55)
print("ECUACION DE ONDA: d²u/dt² = d²u/dx²")
print(f"h = {h},  k = {k},  T = {T}")
print("=" * 55)

# Llamar a la funcion
x, t, U = wave1d(f, g, h, k, T)

m = len(x) - 1
n = len(t) - 1

# --------------------------------------------------------
# Imprimir resultados
# --------------------------------------------------------
print(f"\nVector de nodos espaciales x (tamaño {len(x)}):")
print(np.round(x, 6))

print(f"\nVector de nodos temporales t (tamaño {len(t)}):")
print(np.round(t, 6))

print(f"\nMatriz de aproximaciones U (filas=x, columnas=t), tamaño {U.shape}:")
print(np.round(U, 6))

# Solucion exacta en la malla
U_exacta = np.zeros((m + 1, n + 1))
for i in range(m + 1):
    for j in range(n + 1):
        U_exacta[i, j] = u_exacta(x[i], t[j])

print("\nSolucion exacta en la malla:")
print(np.round(U_exacta, 6))

# Error absoluto maximo
error_max = np.max(np.abs(U_exacta - U))
print(f"\nError absoluto maximo: {error_max:.6e}")

# --------------------------------------------------------
# GRAFICA: comparar solucion exacta vs aproximacion en t = T
# --------------------------------------------------------
plt.figure(figsize=(9, 5))

plt.plot(x, U_exacta[:, -1], 'k-',  linewidth=2, label="Exacta")
plt.plot(x, U[:, -1],        'ro--', linewidth=1.5, markersize=6, label="Aproximacion")

plt.xlabel("Eje x")
plt.ylabel("Solucion u")
plt.title(f"Ecuacion de onda — Diferencias finitas\nt = T = {T},  h = {h},  k = {k},  r = k/h = {k/h:.2f}")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()