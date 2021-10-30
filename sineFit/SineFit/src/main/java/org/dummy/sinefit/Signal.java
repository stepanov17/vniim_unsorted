
package org.dummy.sinefit;

import java.util.Random;

public class Signal {

    private final static Random RND = new Random(100500L);

    public static double[] generate(double x[], double A[], double w[], double phi[], double C) {

        int nw = w.length;
        if (A.length != nw) {
            throw new RuntimeException("nw != nA");
        }
        if (phi.length != nw) {
            throw new RuntimeException("nw != nPhi");
        }

        int nx = x.length;
        double s[] = new double[nx];

        for (int i = 0; i < nx; ++i) {

            s[i] = C;

            for (int j = 0; j < nw; ++j) {

                s[i] += A[j] * Math.cos(w[j] * x[i] + phi[j]);
            }
        }

        return s;
    }

    public static double[] quantization(double x[], double h) {

        int nx = x.length;
        double v[] = new double[nx];

        double c = 1. / h;
        for (int i = 0; i < nx; ++i) { v[i] = h * Math.round(c * x[i]); }

        return v;
    }

    public static double[] addNoise(double x[], double sigma) {

        int nx = x.length;
        double v[] = new double[nx];

        for (int i = 0; i < nx; ++i) { v[i] = x[i] + sigma * RND.nextGaussian(); }

        return v;
    }

    public static double maxDiff(double x[], double y[]) {

        int nx = x.length;
        if (y.length != nx) { throw new RuntimeException("nx != ny"); }

        double d = 0.;
        for (int i = 0; i < nx; ++i) { d = Math.max(d, Math.abs(x[i] - y[i])); }

        return d;
    }

    public static double r2(double x[], double y[]) {

        int nx = x.length;
        if (y.length != nx) { throw new RuntimeException("nx != ny"); }

        double r = 0.;

        for (int i = 0; i < nx; ++i) { r += (x[i] - y[i]) * (x[i] - y[i]); }

        return Math.sqrt(r);
    }
}
