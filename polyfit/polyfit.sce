clear;
clc;

/////////////////////
// polyfit routine //
/////////////////////

// https://stackoverflow.com/questions/30254429/how-to-make-a-polynomial-approximation-in-scilab

function [coeffs] = polyfit(x, y, deg)

    A = ones(length(x), deg + 1)
    for i = 1 : deg
        A(:, i + 1) = x(:) .^ i;
    end
    coeffs = lsq(A, y(:))';

endfunction


function [v] = polyvals(x, coeffs)

    nc = length(coeffs);
    o = ones(1, nc);
    d = 0 : (nc - 1);

    v = zeros(x);
    for i = 1 : length(x)
        a = (x(i) * o) .^ d;
        v(i) = coeffs * a';
    end

endfunction


function [r] = roundTo(x, m)
    
    d1 = 10.^m;
    d2 = 1. / d1;
    r = d2 * round(d1 * x);
    
endfunction

// example: K_0.95 = K_0.95(n0 / n) - exact coverage factor for AR(1) noise regression bands,
// the 1st n0 samples used for the drift prediction for all the n > n0 data
// x_{k + 1} = a x_k + sqrt{1 - a^2} eps_k, eps_k in N(0, 1). here: a = 0.9, n = 100

data = [
0.15  2.420440196990967
0.16  2.421215534210205
0.17  2.421897411346436
0.18  2.422499179840088
0.19  2.423028469085693
0.20  2.423497676849365
0.21  2.423909664154053
0.22  2.424272060394287
0.23  2.424590587615967
0.24  2.424869060516357
0.25  2.425113201141357
0.26  2.425320148468018
0.27  2.425500392913818
0.28  2.425650119781494
0.29  2.425775051116943
0.30  2.425876140594482
0.31  2.425952434539795
0.32  2.426010608673096
0.33  2.426048755645752
0.34  2.426067829132080
0.35  2.426068782806397
0.36  2.426054477691650
0.37  2.426024913787842
0.38  2.425979137420654
0.39  2.425919055938721
0.40  2.425845623016357
0.41  2.425758838653565
0.42  2.425657749176025
0.43  2.425546169281006
0.44  2.425420284271240
0.45  2.425283908843994
0.46  2.425133228302002
0.47  2.424974918365479
0.48  2.424803256988525
0.49  2.424620151519775
0.50  2.424426555633545
0.51  2.424221515655518
0.52  2.424005985260010
0.53  2.423779010772705
0.54  2.423540592193604
0.55  2.423295497894287
0.56  2.423036098480225
0.57  2.422766208648682
0.58  2.422487735748291
0.59  2.422195911407471
0.60  2.421895503997803
0.61  2.421584606170654
0.62  2.421260356903076
0.63  2.420928478240967
0.64  2.420582294464111
0.65  2.420228481292725
0.66  2.419861316680908
0.67  2.419481754302979
0.68  2.419092655181885
0.69  2.418691158294678
0.70  2.418278217315674
0.71  2.417853832244873
0.72  2.417417049407959
0.73  2.416968822479248
0.74  2.416507244110107
0.75  2.416035175323486
];

x = data(:, 1)';
y = data(:, 2)';

polyDeg = 4;
c = polyfit(x, y, polyDeg);
c = roundTo(c, 5) // polynomial coefficieints

///////////////

plot(x, y, '+')

xx = x(1) : 0.01 : x(length(x));

p = polyvals(xx, c);
plot(xx, p, 'k')

legend("K values", "approximation, polyn. deg. = 4")
title("n = 100, alpha = 0.9")
xlabel("n0 / n")
ylabel("K")

p = polyvals(x, c)

// error
e_abs = max(abs(y - p))
e_rel = 100. * max(abs(y - p) ./ y)