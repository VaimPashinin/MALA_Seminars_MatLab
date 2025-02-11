clc;
clear;
%Размер матрицы
n = 23;
%Элементы над- и поддиагоналей
b = 1;
%Элементы диагонали
d = 2;

e = ones(n, 1);
A = full(spdiags([b*e d*e b*e], -1:1, n, n));
lamb_comp = eig(A);
i = 1:n;
lamb_real = d - 2*abs(b)*cos(i*pi/(n + 1));
delta = abs(lamb_real' - lamb_comp);
%Таблица с результатами вычислений
result = table(lamb_comp, lamb_real', delta);
max_delta = max(delta);