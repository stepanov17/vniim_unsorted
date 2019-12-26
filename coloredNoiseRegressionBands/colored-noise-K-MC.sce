// cov. factor calculation example (1^f^alpha noise, alpha = 0.5 .. 1.2)
// https://doi.org/10.32446/0368-1025it.2019-5-21-24
//
// N = 1500

clear;
clc;

stacksize('max');

rand("seed", 1234512345);

RDATA = fscanfMat("ACFDATA.txt");

// norming factor to make the noise stdev at point equal to 1
KDATA = [
 0.50  3.7839301
 0.55  4.2727770
 0.60  4.8111994
 0.65  5.4011682
 0.70  6.0440809
 0.75  6.7406556
 0.80  7.4908460
 0.85  8.2937891
 0.90  9.1477945
 0.95  10.046439
 1.00  10.998359
 1.05  11.986197
 1.10  13.014996
 1.15  14.067397
 1.20  15.160411
];

i0 = 11; // alpha = 1.

R = RDATA(i0, :); // acf

kS = KDATA(i0, 2);

alpha = R(1);

R = R(2 : length(R));


N  = length(R); // 1500
N0 = 500;

W = zeros(N0, N0);
for i = 1 : N0
    l = 1;
    for j = i : N0
        r = R(l);
        W(i, j) = r;
        W(j, i) = r;
        l = l + 1;
    end
end

WInv  = W^(-1);

T0 = 0 : (N0 - 1);
V = [ones(N0, 1), T0'];
Theta = (V' * WInv * V)^(-1);

T = 0 : (N - 1);

U = zeros(T);
for i = 1 : length(T)
    U(i) = sqrt(Theta(1, 1) + 2 * T(i) * Theta(1, 2) + T(i)^2 * Theta(2, 2));
end

M = Theta * V' * WInv;

// the forecast is obtained by means of generalized LSM
function [forecast] = getForecast(vals)

    s = vals(1 : N0);
    C = M * s';
    forecast = C(1) + C(2) * T;

endfunction


function [res] = inBorders(v, L, H)

    n = length(v);
    for i = 1 : n
        if ((v(i) < L(i)) | (v(i) > H(i))) then
            res = 0;
            return;
        end
    end
    res = 1;
endfunction



// 1 / f^alpha noise generator
function [cn] = coloredNoise(alpha, n)

    x = rand(1, n, "normal");
    X = fft(x);
    n2 = n / 2 + 1;
    freqs = 1 : n2;
    X = X(freqs);
    X(1) = 0.;
    X = X ./ (freqs .^ (0.5 * alpha));
    Xc = conj(X(n2 - 1 : -1 : 2));
    y = real(ifft([X Xc]));
    cn = y(1 : n);
    cn = kS * cn; // norming to obtain a unity stdev at point

endfunction


Z = zeros(1, N);


// linear approximation 
function res = getArg(x1, y1, x2, y2)
    k = (y1 - y2) / (x1 - x2);
    b = y1 - k * x1;
    res = (0.95 - b) / k;
endfunction

////////////////////////////////////////////////////////////////////////////////
// Monte Carlo to estimate K for probablity level of 0.95
function [k] = getK(k1, k2, nSim)

    printfFreq = 10000;

    chkS = zeros(1, nSim); // check if stdev at pt = 1

    p1 = 0;
    p2 = 0;

    UExt1 = k1 * U;
    UExt2 = k2 * U;

    for sim = 1 : nSim

        p = coloredNoise(alpha, N);
        fc  = getForecast(p);

        L1 = fc - UExt1;
        H1 = fc + UExt1;

        dp1 = inBorders(Z, L1, H1);
        p1 = p1 + dp1;

        if (dp1 < 0.99) then
            L2 = fc - UExt2;
            H2 = fc + UExt2;
            p2 = p2 + inBorders(Z, L2, H2);
        else
            p2 = p2 + 1; // k1 < k2, so dp1 = 1 => dp2 = 1
        end

        chkS(sim) = p(N);

        if (modulo(sim, printfFreq) == 0) then
            printf("\t%d trials\n", sim);
        end
    end

    kP = 1. / nSim;
    p1 = kP * p1;
    p2 = kP * p2;

    printf("\n>> check: s = %.3f\n", stdev(chkS));
    printf(">> k1 = %.4f, p1 = %.4f\n", k1, p1);
    printf(">> k2 = %.4f, p2 = %.4f\n", k2, p2);

    k = getArg(k1, p1, k2, p2);

endfunction

// iterative approximation
k = getK(1.8, 2.8, 2.e4)
k = getK(k - 0.20, k + 0.20, 2.e4)
k = getK(k - 0.05, k + 0.05, 2.e4)
k = getK(k - 0.05, k + 0.05, 2.e5)
k = getK(k - 0.02, k + 0.02, 2.e5)
k = getK(k - 0.01, k + 0.01, 2.e5)
k = getK(k - 0.005, k + 0.005, 4.e5)
//k = getK(k - 0.0025, k + 0.0025, 5.e5)
