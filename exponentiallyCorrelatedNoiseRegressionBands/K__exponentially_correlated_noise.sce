// see https://doi.org/10.32446/0368-1025it.2019-5-21-24
// cov. factor calculation example (exponentially correlated noise)

clear;
clc;

////////////////////////////////////////////////////////////////////////////////

function [res] = whiteNoise(T)
    res = rand(T, "normal");
endfunction

// exponentially correlated noise
function [res] = expNoise(T, c)
        w = whiteNoise(T);
        nT = length(T);
        res = zeros(T);
        res(1) = w(1);
        for i = 2 : nT
            res(i) = c * res(i - 1) + sqrt(1 - c^2) * w(i);
        end
endfunction

// ACF
function [res] = getR(n, c)
    res = ones(1, n);
    for i = 2 : n
        res(i) = res(i - 1) * c;
    end
endfunction

////////////////////////////////////////////////////////////////////////////////

N  = 50;
N0 = 20;

// exp. correlated noise params
ca = 0.08;
ce = exp(-ca);

////////////////////////////

T = 1 : N;

hT = 0.01;
TT = 1 : hT : N;

R = getR(N, ce);

// for white noise
//R = zeros(1, N);
//R(1) = 1.;

Zeros = zeros(1, N);
ZZ = zeros(1, length(TT));

///////////////////////////////////////////////////////////////////////////////


function [res] = getW(R, n)

    W = zeros(n, n);
    for i = 1 : n
        l = 1;
        for j = i : n
            W(i, j) = R(T(l));
            W(j, i) = R(T(l));
            l = l + 1;
        end
    end
    res = W;
endfunction


WPart = getW(R, N0);
WAll  = getW(R, N);

WPartInv = WPart^(-1);
WAllInv  = WAll^(-1);

VPart = [ones(N0, 1), T(1 : N0)'];
VAll  = [ones(N, 1), T(1 : N)'];

ThetaPart = (VPart' * WPartInv * VPart)^(-1);
ThetaAll  = (VAll' * WAllInv * VAll)^(-1);

U = zeros(1, length(TT));
for i = 1 : length(TT)
    U(i) = sqrt(ThetaPart(1, 1) + 2 * TT(i) * ThetaPart(1, 2) + TT(i)^2 * ThetaPart(2, 2));
end


Mfc = ThetaPart * VPart' * WPartInv;
Mest = ThetaAll * VAll' * WAllInv;

function [forecast] = getForecast(vals)

    s = vals(1 : N0);
    C = Mfc * s';
    forecast = C(1) + C(2) * TT;

endfunction

function [est] = getDriftEstimation(vals)

    C = Mest * vals';
    est = C(1) + C(2) * TT;

endfunction


////////////////////////////////////////////////////////////////////////////////

// TODO: faster check!
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


printfFreq = 1000;

function [P] = getP(kU, nSim)

    P = 0;
    UExt = kU * U;

    for sim = 1 : nSim

        p = expNoise(T, ce);
        //p = whiteNoise(T);

        fc  = getForecast(p);

        //est = getDriftEstimation(p);

        L = fc - UExt;
        H = fc + UExt;

        P = P + inBorders(ZZ, L, H);
        //P = P + inBorders(est, L, H);

        if (modulo(sim, printfFreq) == 0) then
            printf("\t%.3f, %d, %.5f\n", kU, sim, P / sim);
        end
    end

    kP = 1. / nSim;
    P = kP * P;

endfunction


P = [];

////////////////////////////////////////////////////////////////////////////////
// Monte Carlo to estimate K for probablity level of 0.95

nSim = 1e6; // number of trials
K = [2.404 2.405];

for i = 1 : length(K)

    p = getP(K(i), nSim);
    P = [P, p];

    printf("p = %f\n", p);
end;

[K' P']

function res = getArg(x1, y1, x2, y2)
    k = (y1 - y2) / (x1 - x2);
    b = y1 - k * x1;
    res = (0.95 - b) / k;
endfunction

k = getArg(K(1), P(1), K(2), P(2))
