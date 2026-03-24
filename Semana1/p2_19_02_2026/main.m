m = 6;

% Construir vectores (i va de 0 a m-1 según el enunciado)
i = (0:m-1);

a = 2*(i + 1);                    % largo m
b = (i(1:m-1) + 1) / 3;          % largo m-1  (i=0,...,m-2)
c = i(1:m-1) / 3;                 % largo m-1
d = (i(1:m-2) + 2) / 4;          % largo m-2  (i=0,...,m-3)
e = i(1:m-2) / 4;                 % largo m-2

% Construir vector h
h = (2 * i)';

% Construir la matriz y resolver el sistema
A = pentadiagonal(m, a, b, c, d, e);
x = A \ h;

% Calcular e imprimir el errla "noma euclidiana"
error_val = norm(A*x - h, 2);
fprintf("Error ||Ax - h||_2 = %.6e\n", error_val);
