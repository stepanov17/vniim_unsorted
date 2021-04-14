clear;
clc;

rand('seed', 100500);

// signal form:
// y = C + sum_{i = 1..nw} (  A(i) * cos( w(i) * x + phi(i) )  )
//
// nw = 1: P. Handel, Properties of the IEEE-STD-1057 four-parameter sine wave fit algorithm,
// IEEE Transactions on Instrumentation and Measurement, v. 46, 6, 2000


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
            D(:, 2 * i - 1) = COS;
            D(:, 2 * i    ) = SIN;
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

function [y] = signal(x, A, w, phi, C)

    nx = length(x);
    y = C * ones(1, nx);

    nw = length(w);
    for i = 1 : nw
        y = y + A(i) * cos(w(i) * x + phi(i));
    end

endfunction


function [res] = quantization(y, h)
    res = h * round((1. / h) * y);
endfunction


function [res] = addNoise(y, h)
    res = y + h * rand(1, length(y), "normal");
endfunction

////////////////////////////////////////////////////////////////////////////////
// example

A0 = [1., 0.05, 0.03, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01];
nw = length(A0);
w0 = (1 : nw) * %pi;
phi0 = %pi * ones(1, nw);
C0 = 0.;


x_max = 10.;
nx = 256;
h = x_max / (nx - 1);

x = 0 : h : x_max;

y = signal(x, A0, w0, phi0, C0);

// add noise or/and quantization
y = addNoise(y, 1.e-2);
//y = quantization(y, 5.e-2);
plot(x, y);

////////// 3 parameters per harmonic //////////
[A3, phi3, C3] = fit3params(x, w0, y)
est_3par = signal(x, A3, w0, phi3, C3);
plot(x, est_3par, 'k');

////////// 4 parameters per harmonic //////////

// test the algorithm stability: some deviation for the w0 vector
dw0 = 0.01 * max(w0) * ones(w0);

[A4, w4, phi4, C4, nIt] = fit4params(x, w0 + dw0, y)
est_4par = signal(x, A4, w4, phi4, C4);
plot(x, est_4par, 'r');

plot([0., x_max + 3.], [0., 0.], 'k--');
legend("signal", "3 pph approx.", "4 pph approx.")

// compare (e.g., A, w)

printf("\n\nA: initial, 3-par. appr., 4-par. appr\n");
for i = 1 : nw
    printf("%.5f  %.5f  %.5f\n", A0(i), A3(i), A4(i));
end

printf("\n\nw: initial, 4-par. appr\n");
for i = 1 : nw
    printf("%.5f  %.5f\n", w0(i), w4(i));
end

printf("\n\nphi: initial, 3-par. appr., 4-par. appr\n");
for i = 1 : nw
    printf("%.5f  %.5f  %.5f\n", phi0(i), phi3(i), phi4(i));
end
