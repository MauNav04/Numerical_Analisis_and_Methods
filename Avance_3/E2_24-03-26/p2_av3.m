function p2_av3
  a = [2,1,1;4,-6,0;-2,7,2];
  b = [5;-2;9];

  matAumentada = [a,b]

  m = rows(a)

  [U, x] = eliminacionGaussina (matAumentada, m);

endfunction

function [U, x] = eliminacionGaussina (matAumentada, m)
  for k = 1:m-1

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
