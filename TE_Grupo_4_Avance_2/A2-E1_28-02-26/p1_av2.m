function p1_av2()
  clc;

  % Datos del problema
  a = 0;
  b = 2;
  tol = 1e-10;
  maxit = 10000;
  n = 10;
  repeticiones = 100;

  %----------------------------------
  % Ejecucion unica de cada metodo
  %----------------------------------
  [x_bis, err_bis, it_bis, t_bis] = biseccion(a, b, tol, maxit);
  [x_mod, err_mod, it_mod, t_mod] = biseccion_modificada(a, b, n, tol, maxit);

  %----------------------------------
  % Repeticion 100 veces para tiempo promedio
  %----------------------------------
  tiempos_bis = zeros(repeticiones, 1);
  tiempos_mod = zeros(repeticiones, 1);

  for i = 1:repeticiones
    [~, ~, ~, tiempos_bis(i)] = biseccion(a, b, tol, maxit);
    [~, ~, ~, tiempos_mod(i)] = biseccion_modificada(a, b, n, tol, maxit);
  end

  prom_bis = mean(tiempos_bis);
  prom_mod = mean(tiempos_mod);

  %----------------------------------
  % Mostrar resultados
  %----------------------------------
  fprintf('================ RESULTADOS ================\n\n');

  fprintf('Metodo de biseccion clasico:\n');
  fprintf('  x_k               = %.15f\n', x_bis);
  fprintf('  |f(x_k)|          = %.15e\n', err_bis);
  fprintf('  iteraciones       = %d\n', it_bis);
  fprintf('  tiempo (1 corrida)= %.15e s\n', t_bis);
  fprintf('  tiempo promedio   = %.15e s\n\n', prom_bis);

  fprintf('Metodo modificado de la biseccion:\n');
  fprintf('  x_k               = %.15f\n', x_mod);
  fprintf('  |f(x_k)|          = %.15e\n', err_mod);
  fprintf('  iteraciones       = %d\n', it_mod);
  fprintf('  tiempo (1 corrida)= %.15e s\n', t_mod);
  fprintf('  tiempo promedio   = %.15e s\n\n', prom_mod);
 end

%=========================================================
% Funcion del problema
%=========================================================
function y = f(x)
  y = sin(x).^2 - x.^2 + 1;
end

%=========================================================
% Metodo de biseccion clasico
% Devuelve: x, error, iteraciones, tiempo
%=========================================================
function [x, err, iter, tiempo] = biseccion(a, b, tol, maxit)
  t0 = tic;

  fa = f(a);
  fb = f(b);

  if a >= b
    error('Se requiere que a < b.');
  end

  if abs(fa) < tol
    x = a;
    err = abs(fa);
    iter = 0;
    tiempo = toc(t0);
    return;
  end

  if abs(fb) < tol
    x = b;
    err = abs(fb);
    iter = 0;
    tiempo = toc(t0);
    return;
  end

  if fa * fb > 0
    error('No hay cambio de signo en [a,b]. En éste rango no existen raíces de la función.');
  end

  for iter = 1:maxit
    x = (a + b) / 2;
    fx = f(x);
    err = abs(fx);

    if err < tol
      tiempo = toc(t0);
      return;
    end

    if fa * fx < 0
      b = x;
      fb = fx;
    elseif fa * fx > 0
      a = x;
      fa = fx;
    else
      err = 0;
      tiempo = toc(t0);
      return;
    end
  end

  tiempo = toc(t0);
end

%=========================================================
% Metodo modificado de la biseccion
% Divide [a,b] en n subintervalos por iteracion
% Devuelve: x, error, iteraciones, tiempo
%=========================================================
function [x, err, iter, tiempo] = biseccion_modificada(a, b, n, tol, maxit)
  t0 = tic;

  if a >= b
    error('Se requiere que a < b.');
  end

  if n < 2
    error('Se requiere n >= 2.');
  end

  fa = f(a);
  fb = f(b);

  if abs(fa) < tol
    x = a;
    err = abs(fa);
    iter = 0;
    tiempo = toc(t0);
    return;
  end

  if abs(fb) < tol
    x = b;
    err = abs(fb);
    iter = 0;
    tiempo = toc(t0);
    return;
  end

  x = NaN;
  err = NaN;

  for iter = 1:maxit
    h = (b - a) / n;
    encontrado = false;

    for k = 0:(n-1)
      ck = a + k*h;
      ck1 = a + (k+1)*h;

      fck = f(ck);
      fck1 = f(ck1);

      if abs(fck) < tol
        x = ck;
        err = abs(fck);
        tiempo = toc(t0);
        return;
      end

      if abs(fck1) < tol
        x = ck1;
        err = abs(fck1);
        tiempo = toc(t0);
        return;
      end

      if fck * fck1 < 0
        a = ck;
        b = ck1;
        x = (a + b) / 2;
        err = abs(f(x));
        encontrado = true;
        break;
      end
    end

    if ~encontrado
      error('No se encontro subintervalo con cambio de signo.');
    end

    if err < tol
      tiempo = toc(t0);
      return;
    end
  end

  tiempo = toc(t0);
end
