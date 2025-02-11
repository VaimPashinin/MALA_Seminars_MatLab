clc;
clear;
N = 17;
A = randn(N);

% Вычисляем собственные значения
eigenvalues = eig(A);
step = 0.1;
X = min(real(eigenvalues)) - 1:step:max(real(eigenvalues)) + 1;
Y = min(imag(eigenvalues)) - 1:step:max(imag(eigenvalues)) + 1;

% Инициализируем массив для хранения значений функции f(lambda)
f_values = zeros(length(X),length(Y));

% Вычисляем минимальное сингулярное число для каждого значения lambda
for i = 1:length(X)
    for j = 1:length(Y)
        f_values(i, j) = sigma_min(A - (complex(X(i), Y(j))) * eye(N));
    end
end

% Построение поверхности
figure;
surf(X, Y, f_values');
title('Поверхность функции f(\lambda) = log_{10}(\sigma_{min}(A - \lambda I))');
xlabel('Re(\lambda)');
ylabel('Im(\lambda)');
grid on;

% Строим спектральный портрет
figure;
contour(X, Y, f_values', 'ShowText', 'on');
hold on; % Удерживаем текущий график для добавления собственных значений
plot(real(eigenvalues), imag(eigenvalues), 'ro', 'MarkerSize', 8);
hold off; % Освобождаем график
xlabel('Re(\lambda)');
ylabel('Im(\lambda)');
axis equal; % Устанавливаем равные масштабы по осям
grid on; % Включаем сетку

function min_singular_value = sigma_min(X)
    min_singular_value = log10(min(svd(X))); % Используем SVD для нахождения сингулярных чисел
end