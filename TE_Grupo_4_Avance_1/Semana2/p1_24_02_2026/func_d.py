import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import brentq

# El programa genera un gráfica la cual detiene el progreso del programa.
# Una vez que se cierra la gráfica, el programa continúa y encuentra la
# raíz de la función `F(d)` utilizando el método de Brent.
def func_d(d):
    alpha = 4
    sigmadb = 4
    tlambda = 1
    r = 10
    x1 = 7
    x2 = 6
    
    S = math.pi * r**2
    k = (10 * alpha)/np.log(10)
    gd = (2*S/math.pi) * math.acos(d/(2*r)) - d * math.sqrt(r**2 - d**2/4)
    sigmaR2 = sigmadb**2 / (10*alpha)**2
    sigmaC2 = (gd**2/(2*tlambda*k**2)) * (1/gd + 1/S)
    
    f_d = np.log10(x1/d)/(sigmaR2*np.log(10)) + d*(x2-d)/sigmaC2
    
    return f_d

# --- Vectores ---
d_vals = np.linspace(0.1, 19, 500)
y_vals = np.array([func_d(d) for d in d_vals])  # ✅ notación del profesor

# --- Gráfica ---
plt.figure()
plt.plot(d_vals, y_vals, 'b-', linewidth=1.5)
plt.axhline(0, color='r', linestyle='--', linewidth=1.2)
plt.grid(True)
plt.xlabel('d')
plt.ylabel('F(d)')
plt.title('Funcion F(d)')
plt.show()

# --- Encontrar la raíz ---
d_star = brentq(func_d, 1, 10)
print(f'Distancia estimada: d* = {d_star:.10f}')

residuo = abs(func_d(d_star))
print(f'Residuo |F(d*)| = {residuo:.2e}')

"""
El residuo `|F(d*)|` es una medida de la calidad de la solución:

| Valor             | Interpretación |
|     ---           |       ---      |
| `|F(d*)| = 0`     | Solución exacta (imposible en aritmética de punto flotante) |
| `|F(d*)| ≈ 1e-15` | Solución excelente, cerca de la precisión de máquina |
| `|F(d*)| ≈ 1e-6`  | Solución aceptable según tolerancia |
| `|F(d*)|` grande  | La raíz encontrada no es confiable |
"""