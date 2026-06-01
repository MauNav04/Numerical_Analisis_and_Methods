import math
import numpy as np
import sympy as sp


def _maximo_derivada(f, a, b, orden, muestras=200001):
    """
    Aproxima computacionalmente el maximo de |f^(orden)(x)| en [a,b].

    Estrategia:
    1. Se deriva simbolicamente con SymPy.
    2. Se convierte la derivada a una funcion numerica con lambdify.
    3. Se evalua |f^(orden)(x)| en una malla uniforme de puntos del intervalo.
    4. El maximo de esos valores se usa como aproximacion de M.
    """
    X = sp.symbols('x')
    derivada = sp.diff(f, X, orden)
    derivada_num = sp.lambdify(X, derivada, modules=['numpy'])

    puntos = np.linspace(a, b, muestras)
    valores = np.abs(derivada_num(puntos))

    # Se eliminan posibles valores no finitos, por seguridad numerica.
    valores = np.asarray(valores, dtype=float)
    mascara = np.isfinite(valores)
    if not np.any(mascara):
        raise ValueError('No se pudieron obtener valores finitos de la derivada en el intervalo.')

    puntos_validos = puntos[mascara]
    valores_validos = valores[mascara]
    indice_max = int(np.argmax(valores_validos))

    M = float(valores_validos[indice_max])
    punto_max = float(puntos_validos[indice_max])

    return M, punto_max, derivada


def cota_interpolacion(f, a, b, xv, xb):
    """
    Calcula una cota del error de interpolacion |f(xb) - p_n(xb)|.

    Parametros:
        f  : expresion simbolica de SymPy que representa la funcion f(x)
        a  : extremo inferior del intervalo de analisis
        b  : extremo superior del intervalo de analisis
        xv : vector/lista con los nodos de interpolacion [x0, x1, ..., xn]
        xb : punto donde se desea calcular la cota del error

    Retorna:
        ct : cota del error de interpolacion en xb

    Formula utilizada:
        |f(xb) - p_n(xb)| <= M/(n+1)! * |(xb-x0)(xb-x1)...(xb-xn)|

    donde:
        M = max |f^(n+1)(x)| para x en [a,b]
    """
    if a >= b:
        raise ValueError('Debe cumplirse que a < b.')
    if xb < a or xb > b:
        raise ValueError('El punto xb debe estar dentro del intervalo [a,b].')

    xv = np.asarray(xv, dtype=float)
    n = len(xv) - 1
    orden = n + 1

    M, _, _ = _maximo_derivada(f, a, b, orden)

    producto = 1.0
    for xi in xv:
        producto *= (xb - xi)

    ct = (M / math.factorial(orden)) * abs(producto)
    return ct


if __name__ == '__main__':
    X = sp.symbols('x')

    # Datos tomados de la Pregunta 1
    f = sp.log(sp.asin(X)) / sp.log(X)
    a = 0.1
    b = 0.8
    xv = np.arange(0.1, 0.9, 0.1)  # [0.1, 0.2, ..., 0.8]
    xb = 0.55

    n = len(xv) - 1
    orden = n + 1

    M, punto_max, derivada = _maximo_derivada(f, a, b, orden)
    ct = cota_interpolacion(f, a, b, xv, xb)

    producto = np.prod(xb - xv)

    print('--- Pregunta 2: Cota de error del polinomio de interpolacion ---')
    print('Funcion f utilizada:')
    print('f(x) = ln(asin(x)) / ln(x)')
    print()
    print(f'Intervalo de analisis: [{a}, {b}]')
    print('Vector de nodos de interpolacion xv:')
    print(xv)
    print(f'Grado del polinomio interpolador: n = {n}')
    print(f'Orden de la derivada utilizada: n + 1 = {orden}')
    print(f'Punto de evaluacion: xb = {xb}')
    print()
    print('Maximo aproximado de |f^(n+1)(x)| en [a,b]:')
    print(f'M ≈ {M:.12e}')
    print(f'Este maximo se encontro aproximadamente en x = {punto_max:.6f}')
    print()
    print('Producto |(xb-x0)(xb-x1)...(xb-xn)|:')
    print(f'|producto| ≈ {abs(producto):.12e}')
    print()
    print('Cota de error obtenida:')
    print(f'ct ≈ {ct:.12e}')
    print(f'ct ≈ {ct:.12f}')
