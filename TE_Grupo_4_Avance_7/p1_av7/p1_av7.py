# p1_av7.py
# Pregunta 1 - Cota de error para la regla compuesta de Simpson

import math
import numpy as np


def cuarta_derivada_aprox(f, x, delta):
    """
    Aproxima numericamente la cuarta derivada de f en x usando
    una formula centrada de diferencias finitas.
    """
    return (
        f(x - 2 * delta)
        - 4 * f(x - delta)
        + 6 * f(x)
        - 4 * f(x + delta)
        + f(x + 2 * delta)
    ) / (delta ** 4)


def alpha_max_aprox(f, a, b, puntos=2000):
    """
    Aproxima alpha_max = max |f^(4)(x)| en [a,b].
    Retorna una tupla: (alpha_max_aproximado, x_donde_ocurre).
    """
    if a == b:
        return 0.0, a

    a0 = min(a, b)
    b0 = max(a, b)
    longitud = b0 - a0

    # Paso auxiliar para la derivada numerica.
    # Para el ejemplo del enunciado, longitud/100 = 0.005.
    delta = max(longitud / 100.0, 1e-4)

    xs = np.linspace(a0, b0, puntos)
    valores = []

    for x in xs:
        valor = abs(cuarta_derivada_aprox(f, float(x), delta))
        valores.append(valor)

    indice = int(np.argmax(valores))
    return float(valores[indice]), float(xs[indice])


def cota_error_simpson(alpha_max, a, b, n):
    """
    Calcula la cota dada en el enunciado:
    E <= ((b-a) h^4 / 2880) alpha_max, con h = (b-a)/n.
    """
    if n <= 0:
        raise ValueError("n debe ser positivo.")

    longitud = abs(b - a)
    h = longitud / n
    return (longitud * h**4 / 2880.0) * alpha_max


def cota_simpson_puntos(f, a, b, tol):
    """
    Retorna un entero par n tal que la cota de error de Simpson
    compuesta sea menor que la tolerancia tol.

    Entrada:
        f   : funcion a integrar
        a,b : extremos del intervalo
        tol : tolerancia positiva

    Salida:
        n : numero par de subintervalos
    """
    if tol <= 0:
        raise ValueError("La tolerancia debe ser positiva.")

    alpha_max, _ = alpha_max_aprox(f, a, b)

    # Si la cuarta derivada es aproximadamente cero, Simpson integra
    # exactamente polinomios hasta grado 3. Se retorna el menor n par util.
    if alpha_max == 0:
        return 2

    longitud = abs(b - a)

    # De ((b-a)^5 * alpha_max)/(2880*n^4) < tol se obtiene:
    n_teorico = ((longitud**5 * alpha_max) / (2880.0 * tol)) ** 0.25

    n = math.ceil(n_teorico)

    # Simpson compuesto requiere n par.
    if n < 2:
        n = 2
    if n % 2 != 0:
        n += 1

    # Verificacion final estricta por si el redondeo queda justo en la tolerancia.
    while cota_error_simpson(alpha_max, a, b, n) >= tol:
        n += 2

    return n


def simpson_compuesta(f, a, b, n):
    """
    Aproxima la integral de f en [a,b] con Simpson compuesto.
    n debe ser par.
    """
    if n % 2 != 0:
        raise ValueError("Para Simpson compuesto, n debe ser par.")

    h = (b - a) / n
    suma = f(a) + f(b)

    # Terminos con coeficiente 4: indices impares
    for i in range(1, n, 2):
        suma += 4 * f(a + i * h)

    # Terminos con coeficiente 2: indices pares internos
    for i in range(2, n, 2):
        suma += 2 * f(a + i * h)

    return (h / 3.0) * suma


if __name__ == "__main__":
    # Funcion del enunciado: f(x) = e^x (26 - 10x + x^2)
    def f(x):
        return math.exp(x) * (26 - 10 * x + x**2)

    a = 5.0
    b = 5.5
    tol = 1e-8

    n = cota_simpson_puntos(f, a, b, tol)
    h = (b - a) / n
    alpha_max, x_alpha = alpha_max_aprox(f, a, b)
    cota = cota_error_simpson(alpha_max, a, b, n)
    aproximacion = simpson_compuesta(f, a, b, n)

    print("Resultados para la regla compuesta de Simpson")
    print("------------------------------------------------")
    print(f"n obtenido = {n}")
    print(f"h = {h:.15f}")
    print(f"alpha_max aproximado = {alpha_max:.15f}")
    print(f"x donde ocurre alpha_max aprox. = {x_alpha:.15f}")
    print(f"cota de error obtenida = {cota:.15e}")
    print(f"aproximacion de la integral = {aproximacion:.15f}")
