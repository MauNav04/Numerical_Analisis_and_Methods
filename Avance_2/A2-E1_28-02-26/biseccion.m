function x = biseccion(a, b, tol, maxit)
  % Aproxima un cero de f en [a,b] por el metodo de la biseccion
  % Se detiene cuando |f(x)| < tol o cuando se alcanza maxit

  fa = f(a);
  fb = f(b);

  if a >= b
    error('Se requiere que a < b.');
  end

  if fa * fb > 0
    error('No hay cambio de signo en [a,b]. El metodo de biseccion no aplica.');
  end

  for it = 1:maxit
    x = (a + b) / 2;
    fx = f(x);

    % Criterio de parada
    if abs(fx) < tol
      return;
    end

    % Actualizacion del intervalo
    if fa * fx < 0
      b = x;
      fb = fx;
    elseif fa * fx > 0
      a = x;
      fa = fx;
    else
      % fx = 0 exacto
      return;
    end
  end

  warning('Se alcanzo maxit sin cumplir |f(x)| < tol.');
end
