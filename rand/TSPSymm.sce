clear;
clc;

rand("seed", 1234512345);


// van Dorp, Kotz, The Standard Two-Sided Power Distribution and its Properties. The American Statistician, May 2002, Vol. 56, No. 2

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

    // mkae the variance equal to 1
    var = 0.5 / ((p + 1.) * (p + 2.));
    res = (1. / sqrt(var)) * res;

endfunction

tsp = TSP(1.e6, 2.);

histplot(200, tsp)
mean(tsp)
variance(tsp)

