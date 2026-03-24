import math
import time
import numpy as np

# Definición de las funciones y sus derivadas
def f1(x):
    return np.cos(x) - x
def df1(x):
    return -np.sin(x) - 1

def f2(x):
    return np.sin(x)**2 - x**2 + 1
def df2(x):
    return 2*np.sin(x)*np.cos(x) - 2*x

def f3(x):
    return x * np.exp(x**2) - np.sin(x)**2 + 3*np.cos(x) + 5
def df3(x):
    return np.exp(x**2) + 2*x**2*np.exp(x**2) - 2*np.sin(x)*np.cos(x) - 3*np.sin(x)

def f4(x):
    return np.sin(x) + x*np.cos(x)
def df4(x):
    return np.cos(x) + np.cos(x) - x*np.sin(x)

def f5(x):
    return x**2 * np.exp(x**2) - np.sin(x)**2 + x
def df5(x):
    return 2*x*np.exp(x**2) + 2*x**3*np.exp(x**2) - 2*np.sin(x)*np.cos(x) + 1

def f6(x):
    return (x - 1)**3 - 1
def df6(x):
    return 3*(x - 1)**2

def f7(x):
    return (x**2 - 1) / (x**2 + 1) + 1
def df7(x):
    return (2*x*(x**2 + 1) - (x**2 - 1)*2*x) / (x**2 + 1)**2

# Criterio de parada
TOL = 1e-10
MAX_ITER = 10000
def criterio_parada(xk, xk1, fxk): return abs(xk1 - xk) + abs(fxk) < TOL

# Método de Bisección
def biseccion(f, a, b):
    inicio = time.perf_counter()
    if f(a) * f(b) >= 0: return None, None, 0, 0.0, "Intervalo inválido"
    for i in range(MAX_ITER):
        mid = (a + b) / 2.0
        if abs(f(mid)) < TOL: return mid, abs(f(mid)), i, time.perf_counter() - inicio, "OK"
        a, b = (mid, b) if f(a) * f(mid) < 0 else (a, mid)
    return mid, abs(f(mid)), i, time.perf_counter() - inicio, "OK"

# Método de Newton-Raphson
def newton_raphson(f, df, x0):
    inicio = time.perf_counter()
    for i in range(MAX_ITER):
        fxk, dfxk = f(x0), df(x0)
        if abs(dfxk) < 1e-14: 
            print(f"Advertencia: Derivada cercana a cero en x0 = {x0}.")
            return x0, float('nan'), i, time.perf_counter() - inicio, "Derivada cercana a cero"
        x1 = x0 - fxk / dfxk
        if criterio_parada(x0, x1, f(x1)): return x1, abs(f(x1)), i, time.perf_counter() - inicio, "OK"
        x0 = x1
    return x1, float('nan'), i, time.perf_counter() - inicio, "Convergencia no alcanzada"

# Método de la Secante
def secante(f, x0, x1):
    inicio = time.perf_counter()
    for i in range(MAX_ITER):
        fx0, fx1 = f(x0), f(x1)
        if fx1 - fx0 == 0: return None, None, i, time.perf_counter() - inicio, "División por cero"
        x2 = x1 - fx1 * (x1 - x0) / (fx1 - fx0)
        if criterio_parada(x1, x2, f(x2)): return x2, abs(f(x2)), i, time.perf_counter() - inicio, "OK"
        x0, x1 = x1, x2
    return x2, abs(f(x2)), i, time.perf_counter() - inicio, "OK"

# Método de Steffensen
def steffensen(f, x0):
    inicio = time.perf_counter()
    for i in range(MAX_ITER):
        fx0 = f(x0)
        fx0_fx0 = f(x0 + fx0)
        denom = fx0_fx0 - fx0
        if abs(denom) < 1e-14: return x0, abs(fx0), i, time.perf_counter() - inicio, "División por cero"
        x1 = x0 - fx0**2 / denom
        if criterio_parada(x0, x1, f(x1)): return x1, abs(f(x1)), i, time.perf_counter() - inicio, "OK"
        x0 = x1
    return x1, abs(f(x1)), i, time.perf_counter() - inicio, "OK"

# Método de Falsa Posición
def falsa_posicion(f, a, b):
    inicio = time.perf_counter()
    for i in range(MAX_ITER):
        fa, fb = f(a), f(b)
        denom = fb - fa
        if abs(denom) < 1e-14: return None, None, i, time.perf_counter() - inicio, "División por cero"
        xk = a - fa * (b - a) / denom
        fxk = f(xk)
        if criterio_parada(a, xk, fxk): return xk, abs(fxk), i, time.perf_counter() - inicio, "OK"
        a, b = (xk, b) if fa * fxk < 0 else (a, xk)
    return xk, abs(fxk), i, time.perf_counter() - inicio, "OK"

# Método de Müller
def muller(f, x0, x1, x2):
    inicio = time.perf_counter()
    for i in range(MAX_ITER):
        f0 = f(x0)
        f1_val = f(x1)
        f2_val = f(x2)

        h0, h1 = x1 - x0, x2 - x1
        if abs(h0) < 1e-14 or abs(h1) < 1e-14: return None, None, i, time.perf_counter() - inicio, "Puntos demasiado cercanos"
        delta0, delta1 = (f1_val - f0) / h0, (f2_val - f1_val) / h1
        denom_d = x2 - x0
        if abs(denom_d) < 1e-14: return None, None, i, time.perf_counter() - inicio, "División por cero"
        d = (delta1 - delta0) / denom_d
        b_coef = delta1 + h1 * d
        discriminante = b_coef**2 - 4 * f2_val * d
        if discriminante < 0: return None, None, i, time.perf_counter() - inicio, "Raíz compleja"
        sqrt_disc = math.sqrt(discriminante)
        denom_final = b_coef + sqrt_disc if abs(b_coef + sqrt_disc) >= abs(b_coef - sqrt_disc) else b_coef - sqrt_disc
        if abs(denom_final) < 1e-14: return None, None, i, time.perf_counter() - inicio, "División por cero"
        xk = x2 - 2 * f2_val / denom_final
        if criterio_parada(x2, xk, f(xk)): return xk, abs(f(xk)), i, time.perf_counter() - inicio, "OK"
        x0, x1, x2 = x1, x2, xk
    return xk, abs(f(xk)), i, time.perf_counter() - inicio, "OK"

# Intervalos iniciales para cada función
intervalos = {
    "f1": (-1.0, 1.0),
    "f2": (1.0, 2.0),
    "f3": (-2.0, -1.0),
    "f4": (-3.0, -1.5),
    "f5": (-1.0, -0.3),
    "f6": (1.5, 3.0),
    "f7": (-1.0, 0.5),
}

funciones = {
    "f1": f1,
    "f2": f2,
    "f3": f3,
    "f4": f4,
    "f5": f5,
    "f6": f6,
    "f7": f7,
}

# Tabla comparativa
def imprimir_tabla(nombre_func, resultados, a, b):
    print(f"\n{'='*110}")
    print(f"  Función: {nombre_func}    Intervalo base: [{a}, {b}]")
    print(f"{'='*110}")
    encabezado = f"{'Método':<18} {'Valores Iniciales':<28} {'xk':>18} {'Error':>14} {'Iter':>6} {'Tiempo (s)':>12}  Estado"
    print(encabezado)
    print(f"{'-'*110}")

    for metodo, vals_ini, xk, error, iters, tiempo, estado in resultados:
        xk_str = f"{xk:.10f}" if xk is not None else "N/A"
        err_str = f"{error:.4e}" if error is not None else "N/A"
        print(f"{metodo:<18} {vals_ini:<28} {xk_str:>18} {err_str:>14} {iters:>6} {tiempo:>12.6e}  {estado}")

    print(f"{'='*110}")

def ejecutar_todos():
    for nombre, f in funciones.items():
        a, b = intervalos[nombre]
        x0_nr = a
        x0_sec, x1_sec = a, b
        x0_mull, x1_mull, x2_mull = a, 5*(a+b)/6.0, b

        resultados = []

        # Bisección
        xk, err, it, t, estado = biseccion(f, a, b)
        resultados.append(("Bisección", f"a={a}, b={b}", xk, err, it, t, estado))

        # Newton-Raphson
        xk, err, it, t, estado = newton_raphson(f, df1, x0_nr)  # Ajusta df según la función
        resultados.append(("Newton-Raphson", f"x0={x0_nr}", xk, err, it, t, estado))

        # Secante
        xk, err, it, t, estado = secante(f, x0_sec, x1_sec)
        resultados.append(("Secante", f"x0={x0_sec}, x1={x1_sec}", xk, err, it, t, estado))

        # Steffensen
        xk, err, it, t, estado = steffensen(f, x0_nr)
        resultados.append(("Steffensen", f"x0={x0_nr}", xk, err, it, t, estado))

        # Falsa Posición
        xk, err, it, t, estado = falsa_posicion(f, a, b)
        resultados.append(("Falsa Posición", f"a={a}, b={b}", xk, err, it, t, estado))

        # Müller
        xk, err, it, t, estado = muller(f, x0_mull, x1_mull, x2_mull)
        resultados.append(("Müller", f"x0={a:.4f}, x1={x1_mull:.4f}, x2={b}", xk, err, it, t, estado))

        imprimir_tabla(nombre, resultados, a, b)

# MAIN
if __name__ == "__main__":
    print("="*110)
    ejecutar_todos()
    print("\n¡Ejecución completada!\n")