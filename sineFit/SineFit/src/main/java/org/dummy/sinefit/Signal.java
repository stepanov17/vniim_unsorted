
package org.dummy.sinefit;

public class Signal {

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
}
