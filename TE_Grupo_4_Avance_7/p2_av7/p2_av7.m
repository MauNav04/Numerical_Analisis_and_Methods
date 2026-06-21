function p2_av7()
  f = @(x) x^4;
  a = 1;
  b = 3;
  k = 3;
  n = 10;

  I = cuad_gauss(f, a, b, k)

  I = cuad_gauss_comp(f, a, b, k, n)

  % ============================================================================
  % Funcionn para calcular cuadratura gaussiana simple
  % ============================================================================

  function I = cuad_gauss(f, a, b, k)

  % Traemos los valores de x y w previamente computados (se extraen del repositorio de git indicado).
  % Valores se pueden observar en el archivo nodos_pesos.m .
  % Se colocan en archivo aparte meramente por cuestiones de modularización.
  [x, w] = nodos_pesos(k);

  % Se calcula el t con la formula del cambio de variable para un intervalo diferente de [-1, 1]
  t = ((b-a).*x + (b+a)) ./ 2;

  % Se calcula la aproximación de la cuadratura gaussiana en el intervalo dado.
  I = ((b-a)/2) * sum(w .* arrayfun(f, t));

  endfunction

  % ============================================================================
  % Funcionn para calcular cuadratura gaussiana con método compuesto (subintervalos)
  % ============================================================================

  function I = cuad_gauss_comp(f, a, b, k, n)

    % calculamos h en base a la cantidad de subintervalos
    h = (b-a)/n;
    I = 0;

    % Hacemos un ciclo for para trabajar calcular la cuadratura gaussiana simple de subintervalo.
    for j = 1:n
      aj = a + (j-1)*h;
      bj = a + j*h;

      I = I + cuad_gauss(f, aj, bj, k);
    endfor

  endfunction
end
