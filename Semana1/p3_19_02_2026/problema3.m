% =========================================================
% Problema 3 - Octave GUI
% =========================================================
clear; clc; close all;
pkg load symbolic
syms x

% -----------------------------
% Funcion
% -----------------------------
f = (x^3 - 3*x^2 + 3*x - 1)/(x^2 - 2*x);
f = simplify(f);

disp("==============================================");
disp("f(x) = "); pretty(f)

% -----------------------------
% a) Dominio
% -----------------------------
disp("==============================================");
disp("a) Dominio");

[num, den] = numden(f);

disp("Denominador factorizado:");
pretty(factor(den))

poles = solve(den == 0, x);
poles = poles(imag(poles)==0);   % solo reales

disp("Puntos excluidos (den=0):");
disp(poles)
disp("Dominio: R \\ {0, 2}");

% -----------------------------
% b) Intersecciones
% -----------------------------
disp("==============================================");
disp("b) Intersecciones con ejes");

roots_num = solve(num == 0, x);
roots_num = roots_num(imag(roots_num)==0);  % solo reales

disp("Cortes con eje x (num=0):");
disp(roots_num)

disp("Corte con eje y:");
disp("No existe porque x=0 no pertenece al dominio.");

% -----------------------------
% c) Asintotas
% -----------------------------
disp("==============================================");
disp("c) Asintotas");

disp("Asintotas verticales (den=0):");
disp(poles)

% Division polinomial numerica con deconv para asintota oblicua
numP = sym2poly(expand(num));
denP = sym2poly(expand(den));
[qP, rP] = deconv(numP, denP);

q = simplify(poly2sym(qP, x));
r = simplify(poly2sym(rP, x));

disp("Asintota oblicua (cociente q(x)):");
pretty(q)

disp("Residuo r(x):");
pretty(r)

disp("Asintota horizontal:");
disp("No existe (grado(num)=3 > grado(den)=2).");

% -----------------------------
% d) Derivadas (simbolicas)
% -----------------------------
disp("==============================================");
disp("d) Derivadas simbolicas");

fp  = simplify(diff(f, x));
fpp = simplify(diff(fp, x));

disp("f'(x) = ");  pretty(fp)
disp("f''(x) = "); pretty(fpp)

disp("f'(x) factorizada = ");  pretty(factor(fp))
disp("f''(x) factorizada = "); pretty(factor(fpp))

% -----------------------------
% e) Graficas de f, f', f''
% -----------------------------
disp("==============================================");
disp("e) Graficas");

fh   = matlabFunction(f);
fph  = matlabFunction(fp);
fpph = matlabFunction(fpp);

% rangos evitando x=0 y x=2
x1 = linspace(-8, -0.2, 1500);
x2 = linspace( 0.2,  1.8, 1500);
x3 = linspace( 2.2,  8.0, 1500);

figure(1); clf; hold on; grid on;
plot(x1, fh(x1)); plot(x2, fh(x2)); plot(x3, fh(x3));
title("f(x)"); xlabel("x"); ylabel("f(x)");
yl = ylim(); plot([0 0], yl, "--"); plot([2 2], yl, "--");
% asintota oblicua (opcional)
xx = linspace(-8, 8, 400);
plot(xx, xx - 1, "--");
hold off;

figure(2); clf; hold on; grid on;
plot(x1, fph(x1)); plot(x2, fph(x2)); plot(x3, fph(x3));
title("f'(x)"); xlabel("x"); ylabel("f'(x)");
yl = ylim(); plot([0 0], yl, "--"); plot([2 2], yl, "--");
hold off;

figure(3); clf; hold on; grid on;
plot(x1, fpph(x1)); plot(x2, fpph(x2)); plot(x3, fpph(x3));
title("f''(x)"); xlabel("x"); ylabel("f''(x)");
yl = ylim(); plot([0 0], yl, "--"); plot([2 2], yl, "--");
hold off;

% -----------------------------
% f) Creciente / decreciente
% -----------------------------
disp("==============================================");
disp("f) Creciente / decreciente");

crit = solve(fp == 0, x);
crit = crit(imag(crit)==0);  % solo reales

disp("Puntos criticos reales (f'(x)=0):");
disp(crit)

% Resultado (conclusion por analisis de signo de f'(x) factorizada)
disp("Creciente en: (-inf, 1 - sqrt(3)) U (1 + sqrt(3), inf)");
disp("Decreciente en: (1 - sqrt(3), 0) U (0, 2) U (2, 1 + sqrt(3))");
disp("Nota: x=1 es critico pero NO cambia el signo (factor (x-1)^2).");

% -----------------------------
% g) Concavidad
% -----------------------------
disp("==============================================");
disp("g) Concavidad");

fpp_fac = factor(fpp);
disp("f''(x) factorizada:");
pretty(fpp_fac)

disp("Ceros reales de f''(x):");
infl = solve(fpp_fac == 0, x);
infl = infl(imag(infl)==0);  % filtra complejos
disp(infl)

disp("Concava hacia abajo en: (-inf, 0) U (1, 2)");
disp("Concava hacia arriba en: (0, 1) U (2, inf)");
disp("Punto de inflexion: x = 1 (y f(1)=0)");

disp("==============================================");
disp("Fin.");
