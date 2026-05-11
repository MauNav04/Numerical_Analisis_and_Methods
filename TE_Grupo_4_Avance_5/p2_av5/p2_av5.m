% Resolución mediante el algoritmo de Thomas

% ==================================== NOTA ====================================
% Los RESULTADOS DE LA CONSOLA se guardan en un archivo de texto externo llamado
% "resultados_edo2_.txt" ya que son demasiado extensos y en la consola no se
% pueden apreciar correctamente.
% ==============================================================================

function p2_av5()
  diary("resultados_edo2.txt");
  diary on;
  % Definición de las funciones del problema
  p = @(x) -1 ./ x;
  q = @(x) 1 ./ (4 .* x.^2) - 1;
  r = @(x) 0 .* x;  % r(x) = 0 para este problema

  % Condiciones de frontera y dominio
  a = 1; b = 6;
  y0 = 1; yn = 0;

  % Solución exacta
  y_exacta = @(x) sin(6 - x) ./ (sin(5) .* sqrt(x));

  % Valores de h a utilizar
  H = [1, 0.5, 0.2, 0.1, 0.01];

  % Preparar la gráfica
  figure('visible', 'off');
  hold on;

  % Colores y marcadores para cada h
  colores = {'r-o', 'b-s', 'g-d', 'm-^', 'c-v'};

  % Tabla de errores
  fprintf('=============================================================\n');
  fprintf('  RESULTADOS DEL MÉTODO DE DIFERENCIAS FINITAS\n');
  fprintf('=============================================================\n\n');

  errores_max = zeros(1, length(H));

  for k = 1:length(H)
    h = H(k);
    fprintf('-------------------------------------------------------------\n');
    fprintf('  h = %.2f\n', h);
    fprintf('-------------------------------------------------------------\n');

    % Llamada a la función edo2
    [x_vec, y_vec] = edo2(p, q, r, h, a, b, y0, yn);

    % Solución exacta en los nodos
    y_ex = y_exacta(x_vec);

    % Error absoluto máximo
    error_max = max(abs(y_vec - y_ex));
    errores_max(k) = error_max;

    % Imprimir resultados
    fprintf('\n  Vector de nodos x:\n');
    fprintf('  ');
    for i = 1:length(x_vec)
      fprintf('%.4f  ', x_vec(i));
      if mod(i, 8) == 0 && i < length(x_vec)
        fprintf('\n  ');
      end
    end
    fprintf('\n');

    fprintf('\n  Aproximación numérica y:\n');
    fprintf('  ');
    for i = 1:length(y_vec)
      fprintf('%.6f  ', y_vec(i));
      if mod(i, 6) == 0 && i < length(y_vec)
        fprintf('\n  ');
      end
    end
    fprintf('\n');

    fprintf('\n  Solución exacta y_exacta:\n');
    fprintf('  ');
    for i = 1:length(y_ex)
      fprintf('%.6f  ', y_ex(i));
      if mod(i, 6) == 0 && i < length(y_ex)
        fprintf('\n  ');
      end
    end
    fprintf('\n');

    fprintf('\n  Error absoluto máximo: %.10e\n\n', error_max);

    % Graficar la aproximación (solo si h >= 0.1 para no saturar)
    if h >= 0.1
      plot(x_vec, y_vec, colores{k}, 'LineWidth', 1.2, 'MarkerSize', 5, ...
           'DisplayName', sprintf('h=%.2f', h));
    else
      plot(x_vec, y_vec, colores{k}, 'LineWidth', 1.2, 'MarkerSize', 2, ...
           'DisplayName', sprintf('h=%.2f', h));
    end
  end

  % Graficar la solución exacta
  x_fino = linspace(a, b, 500);
  y_fino = y_exacta(x_fino);
  plot(x_fino, y_fino, 'k-', 'LineWidth', 2, 'DisplayName', 'Exacta');

  % Configurar la gráfica
  xlabel('Eje x');
  ylabel('Solucion y');
  title('Método de Diferencias Finitas - Aproximaciones vs Solución Exacta');
  legend('Location', 'best');
  grid on;
  hold off;

  % Guardar la gráfica
  print('-dpng', '-r200', 'grafica_edo2.png');
  fprintf('Gráfica guardada como grafica_edo2.png\n');

  % Tabla resumen de errores
  fprintf('\n=============================================================\n');
  fprintf('  TABLA RESUMEN DE ERRORES MÁXIMOS\n');
  fprintf('=============================================================\n');
  fprintf('  %-12s | %-8s | %-20s\n', 'h', 'n', 'Error máximo absoluto');
  fprintf('  -------------|----------|----------------------\n');
  for k = 1:length(H)
    n = (b - a) / H(k);
    fprintf('  %-12.2f | %-8d | %.10e\n', H(k), n, errores_max(k));
  end
  fprintf('=============================================================\n');
    diary off;
end


function [x_vec, y_vec] = edo2(p, q, r, h, a, b, y0, yn)
  % Si nos dieran una cantidad de muestras habría que primero calcular el h pero en este caso ya nos lo dan

  % Número de subintervalos
  n = round((b - a) / h);

  % Generar los nodos interiores x1, x2, ..., x_{n-1}
  % (n-1 incógnitas)
  x_int = zeros(n-1, 1);

  for j = 1:(n-1)
    x_int(j) = a + j * h;
  end

  % evaluamos x en q y ahora qx guarda un vector
  % con todas las imágenes de x en relación a q
  px = p(x_int);

  % evaluamos x en p y ahora px guarda un vector
  % con todas las imágenes de x en relación a p
  qx = q(x_int);

  % evaluamos x en r y ahora px guarda un vector
  % con todas las imágenes de x en relación a r
  rx = r(x_int);

  m = n - 1;  % tamaño del sistema

  % Ahora que ya tenemos todos los valores de p(x) y q(x) procedemos a realizar el cálculo
  % de las diagonales de A.
  diag_principal = zeros(m, 1);
  diag_superior  = zeros(m-1, 1);
  diag_inferior  = zeros(m-1, 1);

  for j = 1:m
    diag_principal(j) = 2 + h^2 * qx(j);
  end

  for j = 1:(m-1)
    diag_superior(j) = -(1 - (h/2) * px(j));    % coeficiente de y_{j+1} en ecuación j
    diag_inferior(j) = -(1 + (h/2) * px(j+1));   % coeficiente de y_{j-1} en ecuación j+1
  end

  % Formar la matriz tridiagonal A
  A = diag(diag_principal) + diag(diag_superior, 1) + diag(diag_inferior, -1);

  % Construir el vector del lado derecho d
  d = -h^2 * rx;

  %Ahora se procede a calcular el vector d
  %Como nuestra ecuación no tiene una ecuación r(x), vamos a proceder analizamos los valores
  %alpha y beta ya que son los únicos que podría sumar algo al vector d.
  %Como alpha = 1 y beta = 0 el único factor que aporta es alpha para d1.
  d(1) = d(1) + (1 + (h/2) * px(1)) * y0;
  % Ecuación m: c_m * yn se mueve al lado derecho
  d(m) = d(m) + (1 - (h/2) * px(m)) * yn;

  % Resolver el sistema usando mldivide (equivalente al método de Thomas para tridiag)
  y_int = mldivide(A, d);

  % Ensamblar la solución completa
  x_vec = [a; x_int; b]';
  y_vec = [y0; y_int; yn]';
end
