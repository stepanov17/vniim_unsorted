// calculation of coverage factor for regression bands
// noise: AR(1), two-side power noise distribution at point

clear;
clc;

rand("seed", 1234512345);


// van Dorp, Kotz, The Standard Two-Sided Power Distribution and its Properties. The American Statistician, May 2002, Vol. 56, No. 2
// generate a symm. case with zero mean and unity variance
// p = 1 => uniform disrt.
// p = 2 => triangular distr.
function [res] = TSP(n, p)

    c = 0.5;
    res = zeros(1, n);
    U = rand(1, n, "uniform");
    for i = 1 : n
        u = U(i);
        r = 0.;
        if (u < c) then
            r = c * (u / c)^(1. / p);
        else
            r = 1. - (1. - c) * ((1. - u) / (1. - c))^(1. / p);
        end
        res(i) = r - c;
    end

    // make the variance equal to 1
    var = 0.5 / ((p + 1.) * (p + 2.));

    res = (1. / sqrt(var)) * res;

endfunction

// AR(1)
// x_i = alpha * x_{i - 1} + eps_i, eps_i from TSP(p)
function [res] = AR(n, alpha, p)

    // eliminate an influence of the initial state
    // 1.e-15 = alpha^n0  =>  n0 = -15 / log10(alpha)
    n0 = round(-15. / log10(alpha));
    n1 = n0 + n;

    //printf("\nn0 = %d\n", n0);

    res = zeros(1, n1);
    tsp = TSP(n1, p);

    res(1) = tsp(1);

    for i = 2 : n1
       res(i) = alpha * res(i - 1) + tsp(i); 
    end

    res = res(n0 + 1 : n1) * sqrt(1. - alpha^2);

endfunction

function [r] = ACF(n, alpha)
    r = ones(1, n);
    for i = 2 : n
        r(i) = alpha * r(i - 1);
    end
endfunction

////////////////////////////////////////////////////////////////////////////////

N  = 100
N0 = 20


// AR(1) parameters
noiseAlpha = 0.8;
noiseP = 2.;

R = ACF(N, noiseAlpha);

////////////////////////////////////////////////////////////////////////////////

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

        p = AR(N, noiseAlpha, noiseP);
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
            printf("\t%d\n", sim);
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
k = getK(k - 0.0025, k + 0.0025, 4.e5)
