% Grafica de F(d)
d_vals = linspace(0.1, 19, 500);
y_vals = arrayfun(@func, d_vals);

figure;
plot(d_vals, y_vals, 'b-', 'LineWidth', 1.5);
hold on;
plot([0.1, 19], [0, 0], 'r--', 'LineWidth', 1.2);
hold off;
grid on;
xlabel('d');
ylabel('F(d)');
title('Funcion F(d)');

% Encontrar la raiz con fzero
d_star = fzero(@func, [1, 19]);
fprintf('Distancia estimada: d* = %.10f\n', d_star);
