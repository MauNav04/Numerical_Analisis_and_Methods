function p1_av6
  clc; close all;
  diary resultados.txt
  diary on

  % --- Datos del ejercicio ---
  x = 0.1:0.1:0.8;                       % nodos en [0.1, 0.8]
  f = @(t) log(asin(t)) ./ log(t);       % f(x) = log_x(asin(x))
  y = zeros(1, length(x));
  for i = 1:length(x)
    y(i) = f(x(i));
  endfor

  % --- Polinomio interpolador con ddn ---
  [tabla, coef, P] = ddn(x, y);

  % --- Salidas exigidas ---
  printf('\n--- Nodos de interpolacion utilizados ---\n');
  disp(x);
  printf('--- Valores de la funcion evaluados en los nodos ---\n');
  disp(y);
  printf('--- Tabla de diferencias divididas ---\n');
  disp(double(tabla));
  printf('--- Coeficientes del polinomio interpolador ---\n');
  disp(double(coef));
  printf('--- Polinomio interpolador de Newton ---\n');
  disp(vpa(P, 6));
  printf('*** Los resultados se puede apreciar mejor en el archivo "resultados.txt" ***\n');

  % --- Grafica ---
  xx = linspace(0.1, 0.8, 400);
  yf = f(xx);
  Pf = matlabFunction(P);                 % convertir P a funcion evaluable
  yp = Pf(xx);

  figure;
  plot(x, y, 'ko', 'MarkerSize', 8, 'MarkerFaceColor', [0 0.6 0]); hold on;
  plot(xx, yf, 'b-',  'LineWidth', 2);
  plot(xx, yp, 'r--', 'LineWidth', 2);
  hold off;  grid on;
  title('Interpolacion por Diferencias Divididas de Newton');
  xlabel('Eje x');  ylabel('Eje y');
  legend('Datos', 'Funcion exacta', 'Polinomio interpolador');

  diary off
endfunction

function [tabla, coef, P] = ddn(x, y)
  pkg load symbolic
  warning('off', 'OctSymPy:sym:rationalapprox');

% DDN  Interpolacion por Diferencias Divididas de Newton (version SIMBOLICA).
%
%   [tabla, coef, P] = ddn(x, y) construye la tabla de diferencias
%   divididas de Newton y genera el polinomio interpolador.
%
%   Este procedimiento se realizó basándose en el pdf "interpolacion_ddn" que
%   explica el Método de Interpolación por Diferencias Divididas de Newton. La
%   única diferencia es que el grado se calcula a partir de los x e y dados.
%
%   Requiere el paquete simbolico:  pkg load symbolic
%
%   ENTRADAS:
%     x : vector de nodos de interpolacion  [x0; x1; ...; xn]
%     y : vector de valores asociados       [y0; y1; ...; yn],  yi = f(xi)
%
%   SALIDAS:
%     tabla : matriz simbolica (n+1)x(n+1) con la tabla de diferencias
%             divididas. Columna 1 = orden 0 (los valores y), etc.
%     coef  : vector simbolico con los coeficientes de Newton (primera
%             fila de la tabla): [ f[x0], f[x0,x1], ..., f[x0,...,xn] ].
%     P     : EXPRESION SIMBOLICA del polinomio interpolador en la
%             variable simbolica X (objeto 'sym'). Para evaluarla:
%               double(subs(P, X, 0.25))           % en un punto
%               Pf = function_handle(P);           % como function handle

  % --- Pasar los datos a forma simbolica (evita avisos de redondeo) ------
  x = sym(x(:));
  y = sym(y(:));

  n = length(x) - 1;          % grado del polinomio interpolador

  % =====================================================================
  %  PASO 1: tabla de diferencias divididas (calculada simbolicamente)
  %    f[xi,...,xi+k] = ( f[xi+1,...] - f[xi,...] ) / ( x(i+k) - x(i) )
  % =====================================================================
  tabla = sym(zeros(n+1, n+1));
  tabla(:, 1) = y;            % columna 1: diferencias de orden 0 (f[xi]=yi)

  for j = 2:(n+1)                       % j-1 = orden de la diferencia
    for i = 1:(n+2-j)                   % filas validas en esa columna
      tabla(i, j) = ( tabla(i+1, j-1) - tabla(i, j-1) ) ...
                    / ( x(i+j-1) - x(i) );
    end
  end

  % Coeficientes de Newton = primera fila de la tabla
  coef = tabla(1, :);

  % =====================================================================
  %  PASO 2: construir el polinomio SIMBOLICAMENTE
  %    (sigue el Algorithm 1 del PDF: variable simbolica X, termino
  %     acumulativo, suma de coef(k)*termino)
  % =====================================================================
  syms X                       % variable simbolica del polinomio
  polinomio = coef(1);         % primer termino: f[x0]
  termino   = sym(1);          % termino acumulativo (X - x0)(X - x1)...

  for k = 1:n
    termino   = termino * (X - x(k));          % multiplicar por (X - x_{k-1})
    polinomio = polinomio + coef(k+1) * termino;
  end

  P = expand(polinomio);       % expresion simbolica expandida del polinomio

end
