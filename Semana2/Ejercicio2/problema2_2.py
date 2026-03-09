import sympy as sp
from sympy.plotting import plot

# Variable simbólica
x = sp.symbols('x')

# Función
f = sp.exp(x)

# Órdenes solicitados
ordenes = [1, 2, 3, 5, 8]

polinomios = {}

print("\nPolinomios de Taylor de e^x en x = 0\n")

# Polinomios de Taylor
for n in ordenes:
    Pn = sp.series(f, x, 0, n + 1).removeO()
    Pn = sp.simplify(Pn)
    polinomios[n] = Pn

    print(f"P_{n}(x) =")
    print(sp.pretty(Pn))
    print()

# Crear gráfico base con e^x
p = plot(f, (x, -2, 2), show=False, legend=True)
p[0].label = "e^x"

# Agregar polinomios
for n in ordenes:
    Pn = polinomios[n]
    p_aux = plot(Pn, (x, -2, 2), show=False)
    p_aux[0].label = f"P_{n}(x)"
    p.append(p_aux[0])

# Mostrar gráfica
p.title = "Aproximaciones de Taylor de e^x en [-2,2]"
p.xlabel = "x"
p.ylabel = "y"
p.show()

# Cálculo del error
print("\nErrores E_n = integral desde -2 hasta 2 de (e^x - P_n(x))^2 dx\n")

for n in ordenes:
    Pn = polinomios[n]
    En = sp.integrate((f - Pn)**2, (x, -2, 2))
    En = sp.simplify(En)

    print(f"E_{n} =")
    print(sp.pretty(En))
    print(f"Valor aproximado = {sp.N(En)}")
    print()