function p4_av3
  n = 1000;

  A = ones(n, n);
  A(1:n+1:end) = 1001;

  b = ones(n, 1);
  x0 = zeros(n, 1);

  tol = 1e-8;
  max_iter = 1000;

  [~, err, k] = gauss_seidel(A, b, x0, tol, max_iter);

  fprintf('Error final ||Ax-b||_2 = %.12e\n', err);

endfunction


function [x, err, k] = gauss_seidel(A, b, x0, tol, max_iter)
    n = length(b);
    x = x0;

    %Definición de parámetros
    %k es el número de iteraciones
    %i es la fila en la que estamos
    %las sumatorias se hacen por vectoriazación para acelerar el rendimiento.

    for k = 1:max_iter
        x_old = x;   % guardamos la iteración anterior

        %Aquí i es de 1 a n porque tenemos que recorrer las n filas.
        for i = 1:n
            %De cada fila vamos a seleccionar solamente los valores que están por
            %debajo de la diagonal. Esos valores se multiplican por los valores del
            %vector x de la iteración de la FILA anterior para poder ir realizando la sustitución
            %hacia adelante.
            %Cabe denotar que el vector x se actualiza en cada iteración de FILA.
            suma1 = A(i,1:i-1) * x(1:i-1);

            %Nuevamente por cada fila vamos sumando los coeficientes por encima de la diagonal,
            %que corresponden a la parte 𝑈 y lo multiplicamos por el respectivo
            %resultado en x del vector de valores "previo"
            suma2 = A(i,i+1:n) * x_old(i+1:n);

            %Forma final de la fórmula de Gauss-Seidel implementada.
            %Despejamos entonces el valor de x(i), se actualiza en el vector x,
            %este se usará para despejar el x(i+1) para la siguiente fila en
            %la siguiente iteración de FILA.
            x(i) = (b(i) - suma1 - suma2) / A(i,i);
        end

        %Se evalúa si se llegó al error aceptable
        err = norm(A*x - b, 2);

        if err < tol
            return;
        end
    end

    %Se evalúa si se llegó al máximo de iteraciones
    err = norm(A*x - b, 2);
end
