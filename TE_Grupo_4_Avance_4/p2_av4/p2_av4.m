function p2_av4
  clear; clc;
  clc;
  format long g;
  more off;

  diary("resultados_sturm.txt");
  diary on;

  pkg load symbolic
  syms lambda

  % Datos de la matriz tridiagonal
  d = [-1,-3,3,-2,-1,0,0,-3,-2,0,0,1];
  delta = [4,1,-2,3,-3,2,-2,-2,1,3,-1];

  % 1. Polinomio característico por recurrencia de Sturm
  [P, pchar] = sturm_tridiag_sym(d, delta, lambda);

  fprintf("\nPolinomio caracteristico obtenido con Sturm:\n");
  disp(expand(pchar));

  % 2. Intervalo de Gershgorin
  [intervalos, intervalo_unico] = gershgorin_tridiag(d, delta);

  fprintf("\nIntervalos de Gershgorin por fila:\n");
  disp(intervalos);

  fprintf("\nIntervalo unico de Gershgorin:\n");
  disp(intervalo_unico);

  % 3. Convertimos el polinomio simbólico en una función numérica
  f = matlabFunction(pchar, 'vars', lambda);

  % 4. Buscamos subintervalos adecuados dentro de Gershgorin
  n = length(d);
  subintervalos = buscar_subintervalos(f, intervalo_unico, n);

  fprintf("\nSubintervalos donde hay cambio de signo:\n");
  disp(subintervalos);

  % 5. Aplicamos falsa posicion en cada subintervalo
  tol = 1e-10;
  maxit = 100;

  valores_propios = zeros(rows(subintervalos), 1);

  for i = 1:rows(subintervalos)
    a = subintervalos(i, 1);
    b = subintervalos(i, 2);

    valores_propios(i) = falsa_posicion(f, a, b, tol, maxit);
  endfor

  valores_propios = sort(valores_propios);

  fprintf("\nAproximaciones de los valores propios usando falsa posicion:\n");
  disp(valores_propios);

  diary off;
  fprintf("\nLa salida completa fue guardada en resultados_sturm.txt\n");

end

function subintervalos = buscar_subintervalos(f, intervalo, cantidad_raices)

  a = intervalo(1);
  b = intervalo(2);

  m = 1000;
  subintervalos = [];

  while rows(subintervalos) < cantidad_raices

    subintervalos = [];
    x = linspace(a, b, m + 1);

    for i = 1:m
      x1 = x(i);
      x2 = x(i + 1);

      y1 = f(x1);
      y2 = f(x2);

      if y1 * y2 < 0
        subintervalos = [subintervalos; x1, x2];
      endif
    endfor

    if rows(subintervalos) < cantidad_raices
      m = 2*m;
    endif

    if m > 100000
      error("No se encontraron suficientes subintervalos. Aumenta el mallado.");
    endif

  endwhile

endfunction

function [intervalos, intervalo_unico] = gershgorin_tridiag(d, delta)

  d = d(:).';
  delta = delta(:).';

  n = length(d);

  if length(delta) != n - 1
    error("delta debe tener longitud n - 1");
  endif

  R = zeros(1, n);

  if n == 1
    R(1) = 0;
  else
    R(1) = abs(delta(1));
    R(n) = abs(delta(n - 1));

    for i = 2:n-1
      R(i) = abs(delta(i - 1)) + abs(delta(i));
    endfor
  endif

  intervalos = [d(:) - R(:), d(:) + R(:)];
  intervalo_unico = [min(intervalos(:,1)), max(intervalos(:,2))];

endfunction

function [P, pchar] = sturm_tridiag_sym(d, delta, lambda)

  d = sym(d(:).');
  delta = sym(delta(:).');

  n = length(d);

  if length(delta) != n - 1
    error("delta debe tener longitud n - 1");
  endif

  P = cell(n + 1, 1);

  P{1} = sym(1);
  P{2} = expand(d(1) - lambda);

  for k = 2:n
    P{k + 1} = expand((d(k) - lambda)*P{k} - delta(k - 1)^2*P{k - 1});
  endfor

  pchar = P{n + 1};

endfunction

function raiz = falsa_posicion(f, a, b, tol, maxit)

  fa = f(a);
  fb = f(b);

  if fa * fb > 0
    error("No hay cambio de signo en el intervalo dado.");
  endif

  for k = 1:maxit

    c = (a*fb - b*fa)/(fb - fa);
    fc = f(c);

    if abs(fc) < tol
      raiz = c;
      return;
    endif

    if fa * fc < 0
      b = c;
      fb = fc;
    else
      a = c;
      fa = fc;
    endif

  endfor

  raiz = c;

endfunction

