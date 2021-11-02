clear;
clc;

function [b0, b1, t, est] = Deming(x, y, delta)

    mx = mean(x);
    my = mean(y);
    n = length(x);
    c = 1. / (n - 1);
    xc = x - mx;
    yc = y - my;
    sxx = c * sum(xc .^ 2);
    sxy = c * xc * yc';
    syy = c * sum(yc .^ 2);

    b1 = (syy - delta * sxx + sqrt((syy - delta * sxx)^2 + 4 * delta * sxy^2)) / (2 * sxy);
    b0 = my - b1 * mx;

    t = x + (b1 / (b1^2 + delta)) * (y - b0 - b1 * x);
    est = b0 + b1 * t;

endfunction

/// example

m = 15;

ux = 0.2;
uy = 0.3;

ex = ux * rand(1, m, "normal");
ey = uy * rand(1, m, "normal");

x = 1 : m;
y = 0.1 + 0.5 * x;

plot(x, y, 'r--')

x = x + ex;
y = y + ey;

plot(x, y, 'r.')

delta = (uy / ux)^2;

[b0, b1, t, est] = Deming(x, y, delta)
plot(t, est, 'k')

// check: inverse Deming
[c0, c1, est, t] = Deming(y, x, 1 / delta)
plot(t, est, 'k.')
