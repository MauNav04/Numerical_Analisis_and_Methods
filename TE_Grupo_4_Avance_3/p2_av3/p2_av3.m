function p2_av3
  coefs = matDiagonal();

  results = resVector();

  matAumentada = [coefs,results];

  m = rows(coefs);

  %eliminacionGaussina retorna entonces U y x que done
  %U es una matriz superior de coeficientes
  %x es el nuevo vector de resultados
  [U, b] = eliminacionGaussina (matAumentada, m);
  x = sustAtras(U,b,m);

  error_final = norm(coefs*x - results, 2);
  disp(error_final);

endfunction

%Función para crear el vector de resultados segun se solicita
function results = resVector()
  m = 1000;
  results = zeros(m,1);
  for i = 1:m
    results(i)= i/10;
  endfor
endfunction

%Función para crear la matriz diagonal segun se solicita
function coefs = matDiagonal()
  vecbase = ones(1000,1);
  vecDiagPrinc = 5*vecbase;
  vecDiagSec = vecbase(2:end);
  coefs = diag(vecDiagPrinc)+diag(vecDiagSec,1)+diag(vecDiagSec,-1);
endfunction

%Para la sustitución hacia atrás manual el vector en el que se almacenarán los
%resultados ya debe de estar precomputado como un vector de ceros.
function x = sustAtras (U,b,m)
  x = zeros(m,1);

  x(m) = b(m) / U(m,m); % El último elemento

  for i = m-1 : -1 : 1
      suma = 0;
      for j = i+1 : m
          suma = suma + U(i,j) * x(j);
      endfor
      x(i) = (b(i) - suma) / U(i,i);
  endfor

endfunction

function [U, b] = eliminacionGaussina (matAumentada, m)
  for k = 1:m-1

    %======================================================================
    %Funcionaliad de pivote para generalizar la eliminacionGaussina
    %======================================================================
    %Estas funciones se añadieron pero luego se notó que no eran necesarias y relentizaban
    %el código. Se dejan comentadas en caso de ser utiles a futuro.

    %Buscamos un pivote adecuado para resolver la matriz
    %columna_relevante = abs(matAumentada(k:m, k));
    %[valor_max, indice_relativo] = max(columna_relevante);
    %obtenemos el índice de la fila real dentro de la matriz
    #{
    p = indice_relativo + k - 1;

    if valor_max == 0
      error("La matriz no tiene solución única (columna de ceros).");
      break
    endif

    if p ~= k
      matAumentada([k, p], :) = matAumentada([p, k], :);
    endif
    #}

    for i = k+1:m
      m_ik = matAumentada(i,k)/matAumentada(k,k);

      matAumentada(i, k:m+1) = matAumentada(i, k:m+1) - m_ik * matAumentada(k, k:m+1);

    endfor

  endfor

  U = matAumentada(:, 1:m);
  b = matAumentada(:,m+1);

endfunction
