clear;
clc;

rand('seed', 100500);


function [A, phi] = convertToPolar(A0, B0)

    A = sqrt(A0^2 + B0^2);

    phi = 0.; // phi in [0, 2 * %pi)
    if (A0 ~= 0.) then
        if (A0 > 0.) then
            if (B0 >= 0.) then
                phi = atan(B0 / A0);
            else
                phi = atan(B0 / A0) + 2. * %pi;
            end
        else
            phi = atan(B0 / A0) + %pi;
        end
    else
        if (B0 > 0.) then
            phi = %pi / 2;
        else
            phi = 3 * %pi / 2;
        end
    end

endfunction


function [A, phi, C] = fit3params(x, w0, y)

    nx = length(x);
    nw = length(w0);

    A   = zeros(1, nw);
    phi = zeros(1, nw);

    ONES = ones(nx, 1);

    D = zeros(nx, 2 * nw + 1);
    for i = 1 : nw
        D(:, 2 * i - 1) = cos(w0(i) * x)';
        D(:, 2 * i    ) = sin(w0(i) * x)';
    end
    D(:, 2 * nw + 1) = ONES;
    est = linsolve(D' * D, -D' * y');

    for i = 1 : nw

        a = est(2 * i - 1);
        b = est(2 * i);
        [A(i), phi(i)] = convertToPolar(a, -b);
    end

    C = est(2 * nw + 1);

endfunction


function [A, w, phi, C, it] = fit4params(x, w0, y)

    nx = length(x);
    nw = length(w0);

    A   = zeros(1, nw);
    w   = w0;
    phi = zeros(1, nw);

    ONES = ones(nx, 1);

    D = zeros(nx, 2 * nw + 1);
    for i = 1 : nw
        D(:, 2 * i - 1) = cos(w(i) * x)';
        D(:, 2 * i    ) = sin(w(i) * x)';
    end
    D(:, 2 * nw + 1) = ONES;
    est = linsolve(D' * D, -D' * y');

    estPrev = [est; zeros(nw, 1)];

    maxIter = 100;
    EPS = 1.e-10;

    for it = 1 : maxIter

        D = zeros(nx, 3 * nw + 1);
        D(:, 2 * nw + 1) = ONES;

        for i = 1 : nw

            COS = cos(w(i) * x)';
            SIN = sin(w(i) * x)';
            D(:, 2 * i - 1) = cos(w(i) * x)';
            D(:, 2 * i    ) = sin(w(i) * x)';
            a = estPrev(2 * i - 1);
            b = estPrev(2 * i);
            D(:, 2 * nw + 1 + i) = -a * (x' .* SIN) + b * (x' .* COS);

        end

        est = linsolve(D' * D, -D' * y');

        eps = norm(est - estPrev);

        w = w + est((2 * nw + 2) : (3 * nw + 1))';

        estPrev = est;
        if (eps <= EPS) then
            break;
        end

    end

    for i = 1 : nw

        a = est(2 * i - 1);
        b = est(2 * i);
        [A(i), phi(i)] = convertToPolar(a, -b);
    end

    C = est(2 * nw + 1);

endfunction

////////////////////////////////////////////////////////////////////////////////

function [res] = quantization(y, h)
    res = h * round((1. / h) * y);
endfunction

function [res] = addNoise(y, h)
    res = y + h * rand(1, length(y), "normal");
endfunction

function [y] = chirp(x, w0, w1)

    nx = length(x);

    //w0 + (nx - 1) * dw = w1 => dw = (w1 - w0) / (nx - 1)
    dw = (w1 - w0) / (nx - 1);

    w = w0;

    y = zeros(1, nx);

    wtmp = 0.;
    for i = 1 : nx
        y(i) = cos(w * x(i));
        wtmp = w;
        w = w + dw;
    end

endfunction


function [y] = signal(x, A, w, phi, C)

    nx = length(x);
    y = C * ones(1, nx);

    nw = length(w);
    for i = 1 : nw
        y = y + A(i) * cos(w(i) * x + phi(i));
    end

endfunction

////////////////////////////////////////////////////////////////////////////////

x_max = 10;

A0 = [1., 0.4, 0.3];
w0 = [2. * %pi, 7. * %pi, 10.];
phi0 = [%pi / 4., %pi / 3., 0.5];
C0 = 0.;

nx = 256;
h = x_max / (nx - 1);

x = 0 : h : x_max;

//dw0 = 0.01 * rand(w0, "normal");
dw0 = 0.01 * ones(w0);

y = signal(x, A0, w0, phi0, C0);
y = addNoise(y, 0.1);
//y = quantization(y, 0.05);

plot(x, y);

[A, phi, C] = fit3params(x, w0, y)
est_3par = signal(x, A, w0, phi, C);
plot(x, est_3par, 'k');

[A, w, phi, C, it] = fit4params(x, w0 + dw0, y)
est_4par = signal(x, A, w, phi, C);
plot(x, est_4par, 'r');
