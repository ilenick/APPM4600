
clear;
clc;

xvalues = [2 4 6 8 10];

a = 0;
b = 40;

h = 0.001;
n = (b-a)/h;

fprintf('x       Trap Error       Quad Error       Trap Evals    Quad Evals\n');

for x = xvalues

    f = @(t) t.^(x-1).*exp(-t);

    % Composite trapezoidal rule
    t = linspace(a,b,n+1);
    y = f(t);

    I_trap = h*(0.5*y(1) + sum(y(2:end-1)) + 0.5*y(end));

    trap_evals = n+1;

    % Adaptive quadrature
    [I_quad, quad_evals] = quad(f,a,b);

    % Exact value
    I_exact = gamma(x);

    % Relative errors
    err_trap = abs(I_trap-I_exact)/I_exact;
    err_quad = abs(I_quad-I_exact)/I_exact;

    fprintf('%-7d %-16.6e %-16.6e %-13d %d\n', ...
        x,err_trap,err_quad,trap_evals,quad_evals);

end
