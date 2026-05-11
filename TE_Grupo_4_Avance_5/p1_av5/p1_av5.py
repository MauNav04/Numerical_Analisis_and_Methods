import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# PROBLEMA DE CAUCHY
# ============================================================

# --- Definicion de f(x, y) y solucion exacta ---
def f(x, y):
    return (x + y) / x

def y_exacta(x):
    return x * (2 + np.log(x / 2))

a  = 2.0    # extremo izquierdo del intervalo
b  = 10.0   # extremo derecho del intervalo
y0 = 4.0    # condicion inicial y(2) = 4
m  = 7      # numero de puntos (m = n+1, con n=6 subintervalos)
n  = m - 1  # numero de subintervalos

# ============================================================
# PASO 1: Calcular h y los nodos x0, x1, ..., xn
# ============================================================
h = (b - a) / (m - 1)

x = np.zeros(m)
for j in range(m):
    x[j] = a + j * h

print("=" * 60)
print("PROBLEMA DE CAUCHY: y' = (x+y)/x,  y(2)=4,  x en [2,10]")
print(f"Tamano de paso: h = (10-2)/(7-1) = {h:.6f}")
print(f"Nodos x: {x}")
print("=" * 60)

# ============================================================
# PASO 2: Calcular y0, y1, ..., yn con cada metodo numerico
# ============================================================
y_real = np.zeros(m)
for j in range(m):
    y_real[j] = y_exacta(x[j])

# ------------------------------------------------------------
# METODO 1: EULER
# ------------------------------------------------------------
y_euler = np.zeros(m)
y_euler[0] = y0

for j in range(n):
    y_euler[j+1] = y_euler[j] + h * f(x[j], y_euler[j])

# ------------------------------------------------------------
# METODO 2: PREDICTOR-CORRECTOR (Heun)
# ------------------------------------------------------------
y_pc = np.zeros(m)
y_pc[0] = y0

for j in range(n):
    y_pred   = y_pc[j] + h * f(x[j], y_pc[j])
    y_pc[j+1] = y_pc[j] + (h/2) * (f(x[j], y_pc[j]) + f(x[j+1], y_pred))

# ------------------------------------------------------------
# METODO 3: RUNGE-KUTTA ORDEN 2 (punto medio)
# ------------------------------------------------------------
y_rk2 = np.zeros(m)
y_rk2[0] = y0

for j in range(n):
    k1 = f(x[j],       y_rk2[j])
    k2 = f(x[j] + h/2, y_rk2[j] + (h/2)*k1)
    y_rk2[j+1] = y_rk2[j] + h * k2

# ------------------------------------------------------------
# METODO 4: RUNGE-KUTTA ORDEN 3
# ------------------------------------------------------------
y_rk3 = np.zeros(m)
y_rk3[0] = y0

for j in range(n):
    k1 = f(x[j],       y_rk3[j])
    k2 = f(x[j] + h/2, y_rk3[j] + (h/2)*k1)
    k3 = f(x[j] + h,   y_rk3[j] - h*k1 + 2*h*k2)
    y_rk3[j+1] = y_rk3[j] + (h/6) * (k1 + 4*k2 + k3)

# ------------------------------------------------------------
# METODO 5: RUNGE-KUTTA ORDEN 4
# ------------------------------------------------------------
y_rk4 = np.zeros(m)
y_rk4[0] = y0

for j in range(n):
    k1 = f(x[j],       y_rk4[j])
    k2 = f(x[j] + h/2, y_rk4[j] + (h/2)*k1)
    k3 = f(x[j] + h/2, y_rk4[j] + (h/2)*k2)
    k4 = f(x[j] + h,   y_rk4[j] + h*k3)
    y_rk4[j+1] = y_rk4[j] + (h/6) * (k1 + 2*k2 + 2*k3 + k4)

# ------------------------------------------------------------
# METODO 6: TAYLOR ORDEN 2
# ------------------------------------------------------------
y_tay2 = np.zeros(m)
y_tay2[0] = y0

for j in range(n):
    d1 = f(x[j], y_tay2[j])     # y'  = (x+y)/x
    d2 = 1.0 / x[j]             # y'' = 1/x
    y_tay2[j+1] = y_tay2[j] + h * d1 + (h**2 / 2) * d2

# ------------------------------------------------------------
# METODO 7: ADAMS-BASHFORTH 2 PASOS
# ------------------------------------------------------------
y_ab2 = np.zeros(m)
y_ab2[0] = y_rk4[0]
y_ab2[1] = y_rk4[1]

for j in range(1, n):
    fj   = f(x[j],   y_ab2[j])
    fj_1 = f(x[j-1], y_ab2[j-1])
    y_ab2[j+1] = y_ab2[j] + (h/2) * (3*fj - fj_1)

# ------------------------------------------------------------
# METODO 8: ADAMS-BASHFORTH 3 PASOS
# ------------------------------------------------------------
y_ab3 = np.zeros(m)
y_ab3[0] = y_rk4[0]
y_ab3[1] = y_rk4[1]
y_ab3[2] = y_rk4[2]

for j in range(2, n):
    fj   = f(x[j],   y_ab3[j])
    fj_1 = f(x[j-1], y_ab3[j-1])
    fj_2 = f(x[j-2], y_ab3[j-2])
    y_ab3[j+1] = y_ab3[j] + (h/12) * (23*fj - 16*fj_1 + 5*fj_2)

# ------------------------------------------------------------
# METODO 9: ADAMS-BASHFORTH 4 PASOS
# ------------------------------------------------------------
y_ab4 = np.zeros(m)
y_ab4[0] = y_rk4[0]
y_ab4[1] = y_rk4[1]
y_ab4[2] = y_rk4[2]
y_ab4[3] = y_rk4[3]

for j in range(3, n):
    fj   = f(x[j],   y_ab4[j])
    fj_1 = f(x[j-1], y_ab4[j-1])
    fj_2 = f(x[j-2], y_ab4[j-2])
    fj_3 = f(x[j-3], y_ab4[j-3])
    y_ab4[j+1] = y_ab4[j] + (h/24) * (55*fj - 59*fj_1 + 37*fj_2 - 9*fj_3)

# ============================================================
# IMPRIMIR RESULTADOS
# ============================================================

nombres  = ["Euler", "Pred-Corrector", "RK2", "RK3", "RK4",
            "Taylor Ord.2", "AB 2 pasos", "AB 3 pasos", "AB 4 pasos"]
metodos  = [y_euler, y_pc, y_rk2, y_rk3, y_rk4,
            y_tay2, y_ab2, y_ab3, y_ab4]

for nombre, y_aprox in zip(nombres, metodos):
    print(f"\n--- {nombre} ---")
    print(f"{'j':>3}  {'x_j':>8}  {'y_j aprox':>12}  {'y exacta':>12}  {'|error|':>12}")
    for j in range(m):
        error = abs(y_real[j] - y_aprox[j])
        print(f"{j:>3}  {x[j]:>8.4f}  {y_aprox[j]:>12.6f}  {y_real[j]:>12.6f}  {error:>12.4e}")

# ============================================================
# GRAFICA: solucion exacta + todos los metodos
# ============================================================

x_fino = np.linspace(a, b, 300)
y_fino = y_exacta(x_fino)

plt.figure(figsize=(12, 7))

plt.plot(x_fino, y_fino, 'k-', linewidth=2.5, label="Solucion exacta")
plt.plot(x, y_euler,  'o--', label="Euler")
plt.plot(x, y_pc,     's--', label="Predictor-Corrector")
plt.plot(x, y_rk2,   '^-',  label="RK2")
plt.plot(x, y_rk3,   'v-',  label="RK3")
plt.plot(x, y_rk4,   'D-',  label="RK4")
plt.plot(x, y_tay2,  'p--', label="Taylor Ord.2")
plt.plot(x, y_ab2,   'h:',  label="Adams-Bashforth 2p")
plt.plot(x, y_ab3,   'x:',  label="Adams-Bashforth 3p")
plt.plot(x, y_ab4,   '*:',  label="Adams-Bashforth 4p")

plt.xlabel("x")
plt.ylabel("y")
plt.title("Solucion numerica del Problema de Cauchy\ny' = (x+y)/x,  y(2)=4,  x en [2,10]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()