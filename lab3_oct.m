clc;
clear;
N = 50
x = -1:2/N:1;
d_x = 0.1
disturbance = -d_x + 2*d_x*rand(size(x));
fd_x = cos(x) + disturbance
M = 4
system = [];
for i = 0:M
  column = (x.^i)';
  system = horzcat(system, column);
end;
p = system\fd_x'
cond(system)
norm(cos(x)' - system * p)

plot(x, cos(x))
hold on;
plot(x, system * p)
hold off;
