// cov. factor calculation example for multivariate normal distribution
// N = 1500

clear;
clc;

stacksize('max');

rand("seed", 1234512345);

RDATA = fscanfMat("D:\other\exp_rel\ACFDATA.txt");

i0 = 11; // alpha = 1.
//i0 = 2; // alpha = 0.5

R = RDATA(i0, :);

//kS = KDATA(i0, 2);

alpha = R(1);

R = R(2 : length(R));


N  = length(R);
N0 = 750;

Wall = zeros(N, N);
for i = 1 : N
    l = 1;
    for j = i : N
        r = R(l);
        Wall(i, j) = r;
        Wall(j, i) = r;
        l = l + 1;
    end
end

Chol = chol(Wall);

W = Wall(1 : N0, 1 : N0);
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


// correlated noise generator
function [s] = correlatedNoise(alpha, n)

    x = rand(1, n, "normal");
    s = (Chol' * x')';

endfunction



Z = zeros(1, N);


function res = getArg(x1, y1, x2, y2)
    k = (y1 - y2) / (x1 - x2);
    b = y1 - k * x1;
    res = (0.95 - b) / k;
endfunction

////////////////////////////////////////////////////////////////////////////////
// Monte Carlo to estimate K for probablity level of 0.95
function [k] = getK(k1, k2, nSim)

    chkS = zeros(1, nSim);

    p1 = 0;
    p2 = 0;

    UExt1 = k1 * U;
    UExt2 = k2 * U;

    for sim = 1 : nSim

        p = correlatedNoise(alpha, N);
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
