function p3_av6
  % Resolucion del sistema de ecuaciones normales - Ley de Boyle

  % 1. Cargar los datos experimentales
  % Esto cargará las variables guardadas en el archivo (asumimos que se llaman V y P)
  load('datos.mat');
  % whos / Este es un comando interesante que muestra los nombres de los datos extraídos.

  % Asegurarnos de que V y P sean vectores columna para evitar errores de dimension
  V = vals_V(:);
  P = vals_P(:);

  % Numero total de datos
  N = length(V); % Esto deberia dar 4901

  % 2. Calcular las sumatorias del sistema matricial
  % Estas ecuaciones se calcularon a mano según la materia de la clase.
  % El procedimeinto se descirbe en la parte escrita.
  sum_1_V2 = sum(1 ./ (V.^2));
  sum_1_V  = sum(1 ./ V);
  sum_P_V  = sum(P ./ V);
  sum_P    = sum(P);

  % 3. Se construye la matriz de coeficientes (M) y el vector de resultados (C)
  M = [sum_1_V2, sum_1_V;
       sum_1_V,  N];

  C = [sum_P_V;
       sum_P];

  % 4. Resolver el sistema M * x = C
  solucion = mldivide(M,C);

  % Extraer las constantes
  a = solucion(1);
  b = solucion(2);

  % 5. Cálculo del error cuadratico medio (ECM)
  % Evaluamos el modelo ajustado en cada V_i y comparamos con P_i
  P_modelo = a ./ V + b;          % presiones predichas por el modelo
  residuos = P - P_modelo;        % diferencia dato - modelo
  ECM = (1/N) * sum(residuos.^2);

  fprintf('Constante a: %f\n', a);
  fprintf('Constante b: %f\n', b);
  fprintf('Error cuadratico medio (ECM): %f\n', ECM);

  % 6. Graficar datos experimentales y modelo ajustado
  % Para la curva del modelo usamos un V ordenado y denso, asi se ve suave
  Vord = linspace(min(V), max(V), 500)';  % vector columna 1..50
  Pajuste = a ./ Vord + b;

  figure;
  plot(V, P, 'b.', 'MarkerSize', 6);       % datos experimentales (puntos)
  hold on;
  plot(Vord, Pajuste, 'r-', 'LineWidth', 2); % modelo ajustado (linea)
  hold off;

  xlabel('Volumen V');
  ylabel('Presion P');
  title('Ajuste por minimos cuadrados - Ley de Boyle');
  legend('Datos experimentales', 'Modelo P = a/V + b');
  grid on;

end
