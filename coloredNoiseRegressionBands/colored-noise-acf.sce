clear;
clc;

stacksize('max');

rand("seed", 1234512345);

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

endfunction

// ACF
function [r] = ACF(alpha, n)

    nn = 1.e6;
    knn = 1. ./ nn;

    s = coloredNoise(alpha, nn + n);
    s1 = s(1 : nn);
    d1 = stdev(s1);
    r = zeros(1, n);
    for i = 1 : n
        s2 = s(i : i + nn - 1);
        r(i) = (s1 * s2') / (d1 * stdev(s2));
    end
    r = knn * r;

endfunction


function [v] = smoothDecrease(x) // x must decrease monotonously

    v = x;
    for i = 2 : length(x)
        v(i) = min(v(i), v(i - 1));
    end

endfunction


nr = 1500;

DATA = []; // alpha -> ACF

for alpha = 0.5 : 0.05 : 1.2

    navg = 50;
    r = zeros(1, nr);
    for i = 1 : navg // some averaging...
        r = r + ACF(alpha, nr);
        printf("avg #%d\n", i);
    end
    r = (1. / navg) * r;
    r(1) = 1.;
    r = smoothDecrease(r);

    DATA = [DATA; [alpha, r]];

    printf("alpha = %f: done\n", alpha)
end


fprintfMat("ACFDATA.txt", DATA, "%7.4f");
