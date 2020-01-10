clear;
clc;

stacksize('max');

rand("seed", 1234512345);

// 1 / f^alpha noise
function [cn] = coloredNoise(n, alpha)

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

endfunction

// noise combination: w * coloredNoise + (1 - w) * whiteNoise + linear drift ( = kDrift * t)
function [s] = noise(n, w, alpha, kDrift)

    whiteNoise = rand(1, n, "normal");
    s = w * coloredNoise(n, alpha) + (1. - w) * whiteNoise;
    k = 1. / stdev(s);
    s = k * s; // unity stdev

    T = 0 : n - 1;
    s = s + kDrift * T;

endfunction



// Allan variance
function [res] = varAllan(X, tau)

    nInt = floor(length(X) / tau); // number of tau-subintervals
    if (nInt < 2) then
        res = -1000;
        return;
    end

    nX = nInt * tau;

    means = zeros(1, nInt);
    for i = 0 : nInt - 1
        start = i * tau + 1;
        XX = X(start : start + tau - 1);
        means(i + 1) = mean(XX);
    end

    v = 0;
    for i = 2 : nInt
        v = v + (means(i) - means(i - 1))^2;
    end
    res = (0.5 * v) / nInt;

endfunction


// noise parameters
n = 1000;
alpha = 1.;
w = 0.8;
kDrift = 0.0025;

T = 0 : n - 1;
maxTau = round(T(length(T) / 4));

tau = 1 : maxTau;
vAllan = zeros(tau); // averaged

nAvg = 1000;

for i = 1 : nAvg
    s = noise(n, w, alpha, kDrift);
    v = zeros(tau);
    for t = 1 : maxTau
        v(t) = varAllan(s, t);
    end
    vAllan = vAllan + v;
    //printf(">> averaging: %d\n", i);
end

vAllan = vAllan ./ nAvg;

//plot(tau, vAllan, 'k')
