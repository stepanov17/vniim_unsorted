package org.vniim.arbands;

import java.util.Random;

public class NoiseGenerator {

    private final static Random RND = new Random(/*123456L*/);

    // p > 0 => TSP
    // p < 0 => normal
    private double[] eps(int n, double p) {

        double e[] = new double[n];

        if (p < 0.) {
            for (int i = 0; i < n; ++i) { e[i] = RND.nextGaussian(); }
            return e;
        }

        double c = 0.5;

        // make the variance equal to 1
        double var = 0.5 / ((p + 1.) * (p + 2.)); 
        double kv = 1. / Math.sqrt(var);

        for (int i = 0; i < n; ++i) {
            double u = RND.nextDouble();
            double v = 0.;
            if (u < c) {
                v = c * Math.pow(u / c, 1. / p);
            } else {
                v = 1. - (1. - c) * Math.pow((1. - u) / (1. - c), 1. / p);
            }
            e[i] = kv * (v - c);
        }

        return e;
    }

    // AR(1)
    // x_i = alpha * x_{i - 1} + eps_i, eps_i from TSP(p) or normal distr. (p < 0)
    // suppose alpha > 0
    public double[] AR1Noise(int n, double a, double p) {

        double noise[] = new double[n];
        int n0 = 200;
        if (a > 0.) {
            // eliminate an influence of the initial state
            // 1.e-10 = alpha^n0  =>  n0 = -10 / log10(alpha)
            n0 = (int) Math.ceil(-10. / Math.log10(a));
        }

        double b = Math.sqrt(1. - a * a);

        double e[] = eps(n + n0, p);

        double x = e[0];
        for (int i = 1; i < n0; ++i) { x = a * x + b * e[i]; }

        for (int i = 0; i < n; ++i) {
            x = a * x + b * e[n0 + i];
            noise[i] = x;
        }

        return noise;
    }
}
