function p2_av7()
  format long

  f = @(x) cos(x).*exp(x);
  a = 2;
  b = 5;
  k = 7;
  n = 20;

  Isimple = cuad_gauss(f, a, b, k)

  Icompuesto = cuad_gauss_comp(f, a, b, k, n)

  [Iiterativa, err, n] = cuad_gauss_iter(f, a, b, k)

  % ============================================================================
  % Funcion para calcular cuadratura gaussiana simple
  % ============================================================================

  function Is = cuad_gauss(f, a, b, k)

  % Traemos los valores de x y w previamente computados (se extraen del repositorio de git indicado).
  % Valores se pueden observar en el archivo nodos_pesos.m .
  % Se colocan en archivo aparte meramente por cuestiones de modularización.
  [x, w] = nodos_pesos(k);

  % Se calcula el t con la formula del cambio de variable para un intervalo diferente de [-1, 1]
  t = ((b-a).*x + (b+a)) ./ 2;

  % Se calcula la aproximación de la cuadratura gaussiana en el intervalo dado.
  Is = ((b-a)/2) * sum(w .* arrayfun(f, t));

  endfunction

  % ============================================================================
  % Funcion para calcular cuadratura gaussiana con método compuesto (subintervalos)
  % ============================================================================

  function Ic = cuad_gauss_comp(f, a, b, k, n)

    % calculamos h en base a la cantidad de subintervalos
    h = (b-a)/n;
    Ic = 0;

    % Hacemos un ciclo for para trabajar calcular la cuadratura gaussiana simple de subintervalo.
    for j = 1:n
      aj = a + (j-1)*h;
      bj = a + j*h;

      Ic = Ic + cuad_gauss(f, aj, bj, k);
    endfor

  endfunction

  % ============================================================================
  % Funcion para calcular cuadratura gaussiana con método iterativo (con tolerancia)
  % ============================================================================

  function [I, err, n] = cuad_gauss_iter(f, a, b, k)

    % Se definen los valores de tolerancia y subitervalos máximos de acuerdo al enunciado
    tol = 1e-8;
    nmax = 100e6;

    n = 1;
    I_ant = cuad_gauss_comp(f, a, b, k, n);

    % Se calcula la cuadratura gaussiana commpuesta con intervalos de 2*n por ciclo hasta
    % llegar a un nivel de tolerancia aceptable o alcanzar el máximo de subintervalos.
    while n < nmax

      n = 2*n;

      if n > nmax
        n = nmax;
      endif

      I = cuad_gauss_comp(f, a, b, k, n);

      err = abs(I - I_ant);

      if err < tol
        return;
      endif

      I_ant = I;

    endwhile

    warning("Se alcanzó el número máximo de subintervalos sin cumplir la tolerancia.");

  endfunction

end
