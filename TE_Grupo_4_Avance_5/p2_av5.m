function p2_av5
  q = @(x) 1 ./ (4 .* x.^2) - 1;
  [x] = edo2(q, 1,6,1,1,0)
end

function [x] = edo2(q,a,b,h,y0,yn)
  %Si nos dieran una cantidad de muestras habría que primero calcular el h pero
  %en este caso ya nos lo dan

  %se procede con el cálculo de la m
  n = (b-a)/h

  %inciamos mediante el calculo de la x y los intervalos
  %como en total el intervalo va de a hasta b, consideramos la condición de parada n equivalente a b
  x = zeros(1,n);

  %calculamos todos los valores para xj
  for i = 1:n
    x(i) = a + i*h;
  endfor

  %evaluamos x en q y ahora q guarda un vector
  %con todas las imágenes de x en relación a q
  qx = q(x)

  %Calulo de la Matriz A

  a = zeros(1,n-1);
  %b = zeros(1,n-1)
  %c = zeros(1,n-1)

  for j = 1:n-1
    a(j) = 2+h*qx(j);
    %b(j) =
    %c(j) =
  endfor


end



##function [] = edo2(p,q,r,h,a,b,y0,yn)
##end
