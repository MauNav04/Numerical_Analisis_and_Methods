# 19/02/26 - Ejercicio 1

function val = func(d)
  alpha   = 4;    lambda = 1;  r = 10;
  sigmadb = 4;    x1 = 7;      x2 = 6;

  S       = pi * r^2;
  k       = (10 * alpha) / log(10);
  gd      = (2*S/pi) * acos(d/(2*r)) - d * sqrt(r^2 - d^2/4);
  sigmaR2 = sigmadb^2 / (10*alpha)^2;
  sigmaC2 = (gd^2 / (2*lambda*k^2)) * (1/gd + 1/S);

  val = log10(x1/d) / (sigmaR2 * log(10)) + (d*(x2-d)) / sigmaC2;
endfunction
