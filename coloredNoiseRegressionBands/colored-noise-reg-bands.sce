clear;
clc;

ACFDATA = fscanfMat("ACFDATA.txt");
[m, n] = size(ACFDATA);

alpha_start = ACFDATA(1, 1);
alpha_end   = ACFDATA(m, 1);
alpha_step  = ACFDATA(2, 1) - alpha_start;
alpha_range = string(alpha_start) + " : " + string(alpha_step) + " : " + string(alpha_end);


// all t points
T = 0 : n - 2;

// uncorrelated white noise ACF
RW = zeros(T);
RW(1) = 1.;


/// plot acf data ///
function plotACFs()

    for i = 1 : m
        acfn = ACFDATA(i, 2 : n);
        plot(T, acfn, 'k');
    end

    title("ACF, 1 / f^alpha noise, alpha = " + alpha_range + " (bottom to top)")

endfunction

////////////
/// GLSM ///
////////////

function [res] = getW(R, k)

    W = zeros(k, k);
    for i = 1 : k
        l = 1;
        for j = i : k
            w = R(T(l) + 1);
            W(i, j) = w;
            W(j, i) = w;
            l = l + 1;
        end
    end
    res = W;

endfunction

// get 2x2 Theta matrix for the GLSM such that the regression band is
// u(t) = sqrt(Theta(1, 1) + 2 * t * Theta(1, 2) + t^2 * Theta(2, 2))
// (for a unity noise dispersion)
function [Theta] = getTheta(R, n0)

    W = getW(R, n0);
    Winv = W^(-1);
    V = [ones(n0, 1), T(1 : n0)'];
    Theta = (V' * Winv * V)^(-1);

endfunction

// range of the data used for the reg. bands estimation: [0 .. n0 - 1]
// noise variance in point is supposed to be equal to 1.
function [u] = getBand(R, n0, t)

    Theta = getTheta(R, n0);

    nt = length(t);
    u = zeros(1, nt);
    for i = 1 : nt
        u(i) = sqrt(Theta(1, 1) + 2. * t(i) * Theta(1, 2) + t(i)^2 * Theta(2, 2));
    end

endfunction

// fund t such that band u(t0) >= u0. the estimation interval for the band is [0 .. n0 - 1]
function [t] = getTByU(R, n0, t0, dt, u0)

    Theta = getTheta(R, n0);
    t = t0;
    while %T
        u = sqrt(Theta(1, 1) + 2. * t * Theta(1, 2) + t^2 * Theta(2, 2));
        if (u >= u0) then
            break;
        end
        t = t + dt;
    end

endfunction

function plotBands(n0, n1)

    dt = 0.1;
    t = 0 : dt : n1 - 1.;

    for i = 1 : m
        R = ACFDATA(i, 2 : n);
        u = getBand(R, n0, t);
        plot(t, u, 'k');
    end

    // uncorrelated white noise case
    u = getBand(RW, n0, t);
    plot(t, u, 'r');
    plot([n0 n0], [0., 1.5], 'k--');

    title("reg. bands, 1 / f^alpha noise, alpha = " + alpha_range + " (bottom to top), red is for white noise. n0 = " + string(n0));

endfunction

// t_K is such a value for t that u(t_K) ~ K * u(n0 - 1)
function [TK] = getTK(K, n0, dt)

    // white noise
    u_tmp = getBand(RW, n0, [n0 - 1]);
    u0 = K * u_tmp(1);
    TK = [getTByU(RW, n0, n0 + dt, dt, u0)];

    for i = 1 : m
        R = ACFDATA(i, 2 : n);
        u_tmp = getBand(R, n0, [n0 - 1]);
        u0 = K * u_tmp(1);
        t = getTByU(R, n0, n0 + dt, dt, u0);
        TK = [TK; t]
    end

endfunction

function [TK] = getTKTable(n0)

    dt = 0.01;

    TK = [-999; ACFDATA(:, 1)];
    for K = [2 3 5]
        TK = [TK, getTK(K, n0, dt)];
    end

endfunction

//plotACFs()
plotBands(200, 500)

getTKTable(100)
getTKTable(200)
getTKTable(300)
getTKTable(500)
