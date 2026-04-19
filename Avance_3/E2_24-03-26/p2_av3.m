function p2_av3
  a = [0,2,1;1,1,2;2,0,1];
  b = [3;2;1];

  matAumentada = [a,b]

  m = rows(a)

  [U, x] = eliminacionGaussina (matAumentada, m);

endfunction

function [U, x] = eliminacionGaussina (matAumentada, m)
  for k = 1:m-1

    %Buscamos un pivote adecuado para resolver la matriz
    columna_relevante = abs(matAumentada(k:m, k));
    [valor_max, indice_relativo] = max(columna_relevante);
    %obtenemos la final real dentro de la matriz
    p = indice_relativo + k - 1;

    if valor_max == 0
      error("La matriz no tiene solución única (columna de ceros).");
      break
    endif

    if p ~= k
      matAumentada([k, p], :) = matAumentada([p, k], :);
    endif

    for i = k+1:m
      m_ik = matAumentada(i,k)/matAumentada(k,k);

      for j = k:m+1
        matAumentada(i,j) = matAumentada(i,j) - m_ik*matAumentada(k,j);
      endfor

    endfor

  endfor

  U = matAumentada(:, 1:m)
  x = matAumentada(:,m+1)

endfunction
