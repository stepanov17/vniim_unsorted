// bayesian refinement of unknown variance

clear;
clc;

// samples: from normal distr. with unknown expected value and variance (sigma);
// exp. may vary from sample to sample, var. is assumed to be the same
SAMPLES = [
 [1.531, 1.479, 1.289, 1.333, 1.357, 1.321, 1.356, 1.442, 1.421, 1.422];
 [1.790, 1.642, 1.702, 1.624, 1.780, 1.755, 1.707, 1.703, 1.684, 1.710];
 [1.399, 1.370, 1.329, 1.274, 1.501, 1.387, 1.440, 1.501, 1.398, 1.404];
];

// coverage probability
//P0 = 0.95;
P0 = 0.99;

sz = size(SAMPLES);
nSamples = sz(1);
n = sz(2);

hs = 5.e-5;
sigma = 0. : hs : 1.;
sigma(1) = 1.e-10; // avoid div by 0


function [i1, I] = upperBnd(x, pdf, p0, i0)

    I = 0.;
    dx = x(2) - x(1);
    i1 = -1;

    f1 = pdf(i0);
    for i = i0 + 1 : length(pdf)
        f2 = pdf(i);
        I = I + 0.5 * (f1 + f2) * dx;
        if (I >= p0) then
            i1 = i;
            break;
        end
        f1 = f2;
    end

    if (I < p0) then
        i1 = -1;
        I = -1;
    end

endfunction

// coverage interval
function [x0, x1, p] = covInterval(x, pdf, p0)

    I0 = [];
    I1 = [];
    CP = [];

    for i0 = 1 : length(x)
        [i1, I] = upperBnd(x, pdf, p0, i0);
        if (i1 == -1) then
            break;
        end
        I0 = [I0, i0];
        I1 = [I1, i1];
        CP = [CP, I];
    end

    DI = I1 - I0;
    minDI = min(DI);

    ind = [];
    for i = 1 : length(CP)
        if (DI(i) == minDI) then
            ind = [ind, i];
        end
    end

    TMP = [I0(ind)', I1(ind)', CP(ind)'];
    [v, j0] = min(CP(ind));

    tmp = TMP(j0, :);
    dx = x(2) - x(1);
    x0 = dx * tmp(1);
    x1 = dx * tmp(2);
    p = tmp(3); // just for check (p0 ~ p)

endfunction

S = 0.;

col = ['b', 'm', 'r', 'k'];

for k = 1 : nSamples

    x = SAMPLES(k, :);
    x = x - mean(x); // centered
    S = S + 0.5 * sum(x .^ 2);
    d = 0.5 * k * (n - 1);
    c = (2. * S^d) / gamma(d);
    e = sqrt(S) * gamma(d - 0.5) / gamma(d); // expectation E(sigma)
    E(k) = e;
    d = k * (n - 1) + 1;
    pdf = c * exp(-S * sigma .^ (-2)) .* sigma .^ (-d);
    f = c * exp(-S * e^(-2)) * e^(-d);

    [x0, x1, p] = covInterval(sigma, pdf, P0); // coverage interval for sigma

    // >>>>> plots >>>>>

    step = 5;

    ind = (x0 / hs) : step : (x1 / hs);
    plot(sigma(ind), pdf(ind), col(k), 'LineWidth', 2);

    f0 = c * exp(-S * x0^(-2)) * x0^(-d);
    plot([x0, x0], f0 * [0.7, 1.3], col(k), 'LineWidth', 2);

    f1 = c * exp(-S * x1^(-2)) * x1^(-d);
    plot([x1, x1], f1 * [0.7, 1.3], col(k), 'LineWidth', 2);

    ind = 1 : step : (x0 / hs - 1);
    plot(sigma(ind), pdf(ind), col(k) + '--', 'LineWidth', 2);

    ind = (x1 / hs + 1) : step : min(1.5 * x1 / hs, length(pdf));
    plot(sigma(ind), pdf(ind), col(k) + '--', 'LineWidth', 2);

    plot(e, f, col(k) + '.');

    // <<<<< plots <<<<<

    printf("#%d: E(sigma) = %.5f, cov. interval (sigma) = [%.5f, %.5f], p0 = %.3f\n", k, e, x0, x1, p);

    //hs * sum(pdf) // check: integral of pdf ~ 1.
end

