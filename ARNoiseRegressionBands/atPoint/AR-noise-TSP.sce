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

    res = res(n0 + 1 : n1) * sqrt(1. - alpha^2); // unity variance for noise

endfunction

nSim = 2e5;
s1 = zeros(1, nSim);
s2 = zeros(1, nSim);
s3 = zeros(1, nSim);

// AR(1) parameters
noiseAlpha = 0.8;
noiseP = 2.;

for sim = 1 : nSim
    x = AR(100, noiseAlpha, noiseP);
    s1(sim) = x(10);
    s2(sim) = x(50);
    s3(sim) = x(100);
    if (modulo(sim, 10000) == 0) then
        printf("\t%d\n", sim);
    end
end

// check that variance is equal to 1
variance(s1)
variance(s2)
variance(s3)

// distr. at point
histplot(200, s3)

x = -3 : 0.1 : 3;
f = zeros(x);
for i = 1 : length(x)
    f(i) = exp(-0.5 * x(i)^2);
end
f = (1. / sqrt(2. * %pi)) * f

plot(x, f, 'r')

