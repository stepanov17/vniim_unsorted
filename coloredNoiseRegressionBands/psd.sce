// check: the noise is indeed 1 / f^alpha noise, i.e. its PSD ~ 1 / f^alpha

clear;
clc;

fact = 1.; // norming factor

alpha = 1.;


// 1 / f^alpha noise generator
function [cn] = coloredNoise(alpha, n)

    x = rand(1, n, "normal");
    X = fft(x);
    n2 = n / 2 + 1;
    freqs = 1 : n2;
    X = X(freqs);
//    X(1) = 0.; // PSD(1) = 0 if uncomment this
    X = X ./ (freqs .^ (0.5 * alpha));
    Xc = conj(X(n2 - 1 : -1 : 2));
    y = real(ifft([X Xc]));
    cn = fact * y(1 : n);

endfunction

function [psd] = PSD(s)

    ns = length(s);
    psd = (1. / ns) * abs(fft(s)).^2; // PSD ~ square of FFT magnitude

endfunction

n = 100

psd = zeros(1, n);

navg = 1.e5;
for i = 1 : navg
    s = coloredNoise(alpha, n);
    psd = psd + PSD(s);
end

psd = (1. / (fact^2 * navg)) * psd;
psdRef = 1. ./ (1 : n).^alpha;

T = 1 : 30;

// compare: PSD and 1 / f^alpha
[T' psdRef(T)' psd(T)']
