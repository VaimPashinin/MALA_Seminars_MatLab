clc;
clear;
% Число точек
M = 7;
% Длина шага
h = 1 / M;
% Число интервалов
N = M - 1;
e = ones(N, 1);
C = diag(e/h^2);
B = full(spdiags([e/h^2 -4*e/h^2 e/h^2], -1:1, N, N));
A = kron(diag(ones(N - 1, 1), -1) + diag(ones(N - 1, 1), 1), C) + kron(eye(N), B);
[V, D] = eig(A);
lambda_comp = diag(D);
lambda_real = [];
for i = 1:N
    for j = (1:N)'
        lr = -4*M^2 * (sin(i*pi/(2*M)).^2 + sin(j*pi/(2*M)).^2);
        lambda_real = [lambda_real;lr];
    end
end
lambda_real = sort(lambda_real);
delta = abs(lambda_real - lambda_comp);
% Таблица с результатами
res1 = table(lambda_comp, lambda_real, delta);
max_delta = max(delta);
% Проверка ортогональности
ort = norm(V * V' - eye(N^2));
% Невязка
disc = norm(A*V - V*D);
