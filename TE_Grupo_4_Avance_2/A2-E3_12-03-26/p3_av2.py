import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

def newton_raphson(x0, f_strings, x_strings, tol=1e-10, iterMax=1000):

    m = len(f_strings)
    assert len(x_strings) == m, "Dimensión inconsistente entre variables y funciones."

    # Construir símbolos y expresiones 
    vars_sym = sp.symbols(x_strings)
    F_sym    = [sp.sympify(fi) for fi in f_strings]

    # Calcular Jacobiano
    J_sym = sp.Matrix([[sp.diff(fi, xj) for xj in vars_sym] for fi in F_sym])

    # Lambdificar F y J 
    F_num = sp.lambdify(vars_sym, sp.Matrix(F_sym), modules='numpy')
    J_num = sp.lambdify(vars_sym, J_sym,            modules='numpy')

    xk     = np.array(x0, dtype=float)
    errores = []
    # Ciclo principal de Newton-Raphson
    for k in range(1, iterMax + 1):
        Fk    = np.array(F_num(*xk), dtype=float).flatten()
        Jk    = np.array(J_num(*xk), dtype=float)
        error = np.linalg.norm(Fk, 2)
        errores.append(error)

        if error < tol:
            break

        # Resolver J * delta = -F 
        try:
            delta = np.linalg.solve(Jk, -Fk)
        except np.linalg.LinAlgError:
            print(f"  [!] Jacobiana singular en iteración {k}. Se detiene.")
            break

        xk = xk + delta

    error_final = np.linalg.norm(np.array(F_num(*xk), dtype=float).flatten(), 2)

    # Gráfica de convergencia 
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.semilogy(range(1, len(errores) + 1), errores, 'o-', color='royalblue',
                markersize=6, linewidth=1.8, label=r'$\|F(x^{(k)})\|_2$')
    ax.axhline(y=tol, color='tomato', linestyle='--', linewidth=1.5,
               label=f'Tolerancia = {tol:.0e}')
    ax.set_xlabel('Iteración $k$',          fontsize=12)
    ax.set_ylabel(r'$\|F(x^{(k)})\|_2$',   fontsize=12)
    ax.set_title('Convergencia de Newton-Raphson', fontsize=13)
    ax.legend(fontsize=11)
    ax.grid(True, which='both', linestyle=':', alpha=0.6)
    plt.tight_layout()

    return xk, k, error_final, errores, fig


# Sistemas de prueba

sistemas = [
    {
        "nombre": "Sistema 1 (2x2 — polinomial)",
        "f_strings": [
            "x1**2 - 2*x1 - x2 + 0.5",
            "x1**2 + 4*x2**2 - 4"
        ],
        "x_strings": ["x1", "x2"],
        "x0": [3.0, 2.0],
    },
    {
        "nombre": "Sistema 2 (2x2 — trigonométrico)",
        "f_strings": [
            "sin(x1) + x2*cos(x1)",
            "x1 - x2"
        ],
        "x_strings": ["x1", "x2"],
        "x0": [1.2, -1.5],
    },
    {
        "nombre": "Sistema 3 (4x4 — cúbico simétrico)",
        "f_strings": [
            "x2*x3 + x4*(x2 + x3)",
            "x1*x3 + x4*(x1 + x3)",
            "x1*x2 + x4*(x1 + x2)",
            "x1*x2 + x1*x3 + x2*x3 - 1"
        ],
        "x_strings": ["x1", "x2", "x3", "x4"],
        "x0": [-1.0, -1.0, -1.0, -1.0],
    },
]


# Ejecución principal

def separador(char='=', n=72):
    print(char * n)

def imprimir_resultado(s, xk, iters, err):
    separador()
    print(f"  {s['nombre']}")
    separador()
    print(f"  x inicial      : {s['x0']}")
    print(f"  Iteraciones    : {iters}")
    print(f"  Error final    : ||F(x^(k))||_2 = {err:.6e}")
    print(f"  Solución aprox.:")
    for xi, name in zip(xk, s['x_strings']):
        print(f"    {name} = {xi:.14f}")
    # Verificación
    vars_sym = sp.symbols(s['x_strings'])
    F_num    = sp.lambdify(vars_sym, sp.Matrix([sp.sympify(fi) for fi in s['f_strings']]), modules='numpy')
    Fval     = np.array(F_num(*xk), dtype=float).flatten()
    print(f"  Verificación F(x*):")
    for i, (fi_str, fv) in enumerate(zip(s['f_strings'], Fval)):
        print(f"    f{i+1}(x*) = {fv:.6e}   [{fi_str} = 0]")
    separador()
    print()


if __name__ == "__main__":
    separador('=')
    print("  NEWTON-RAPHSON PARA SISTEMAS NO LINEALES")
    print("  Criterio de parada: ||F(x^(k))||_2 < tol")
    separador('=')
    print()

    for i, s in enumerate(sistemas, start=1):
        xk, iters, err, errores, fig = newton_raphson(
            x0        = s['x0'],
            f_strings = s['f_strings'],
            x_strings = s['x_strings'],
            tol       = 1e-10,
            iterMax   = 1000
        )
        imprimir_resultado(s, xk, iters, err)

        # Título con nombre del sistema y guardar PNG
        fig.suptitle(s['nombre'], fontsize=11, style='italic')
        nombre_archivo = f"convergencia_sistema{i}.png"
        fig.savefig(nombre_archivo, dpi=150, bbox_inches='tight')
        print(f"  Gráfica guardada: {nombre_archivo}")
        print()

    # Mostrar todas las gráficas al final
    plt.show()