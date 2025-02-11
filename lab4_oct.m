clc;
clear;

M = 51;
N = 91;
a = -16;
b = 16;
c = -16;
d = 16;

h_t = (b-a)/(M-1);
h_s = (d-c)/(N-1);

t = (a:h_t:b)';
s = (c:h_s:d)';

A = zeros(N, M);
for j = 1:N
    for i = 1:M
        if i > 1
            fun_1 = @(x) K_1(x,s(j)) .* (x - t(i - 1)) / h_t;
            A(j, i) = A(j, i) + integral(fun_1, t(i - 1), t(i));
        end
        if i < M
            fun_2 = @(x) K_1(x,s(j)) .* (t(i + 1) - x) / h_t;
            A(j, i) = A(j, i) + integral(fun_2, t(i), t(i + 1));
        end
    end
end
mu = cond(A);

F = f1(s, a, b);
B = zeros(M, M);
for j = 1:M
    for i = 1:M
        if abs(i - j) == 0
            B(i, j) = -2/h_t^2;
        end
        if abs(i - j) == 1
            B(i, j) = 1/h_t^2;
        end
    end
end

tau = 10^-1;
C = [(1 - tau)*A; tau*B];
cond(C)
F_reg = [(1 - tau)*F; zeros(M, 1)];
x_comp = C\F_reg;
x_real = true_x_1(t);

plot(t, x_comp)
hold on;
plot(t, x_real)
hold off;

function true_solution_1 = true_x_1(T)
    true_solution_1 = T.^2;
end

function true_solution_2 = true_x_2(T)
    true_solution_2 = exp(-T);
end

function K_first = K_1(t, s)
    K_first = exp(t - s);
end

function K_second = K_2(t, s)
    K_second = (t - s).^2;
end

function right_first = f1(S, a, b)
        right_first = ((b^2 - 2*b + 2)*exp(b) + (-(a^2) + 2*a - 2)*exp(a))*exp(S);
end

function right_second = f2(S, a, b)
        right_second = exp(-a)*(S.^2 + (-2*a - 2)*S + a^2 + 2*a + 2)...
            - exp(-b)*(S.^2 + (-2*b - 2)*S + b^2 + 2*b + 2);
end
