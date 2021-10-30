
package org.dummy.sinefit;

import java.util.Arrays;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.DecompositionSolver;
import org.apache.commons.math3.linear.QRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;


// P. Handel, Properties of the IEEE-STD-1057 four-parameter sine wave fit algorithm,
// IEEE Transactions on Instrumentation and Measurement, v. 46, 6, 2000

public class SineFit {

    private static final int    MAXITER = 30;
    private static final double EPS = 1.e-10;

    private static final double F_Q = 14400; // Hz

    private double[] toPolar(double x, double y) {

        double res[] = new double[2];
        res[0] = Math.sqrt(x * x + y * y);

        double arg = Math.atan2(y, x);
        if (arg < 0.) { arg += 2. * Math.PI; } // arg in [0, 2 * pi]
        res[1] = arg;

        return res;
    }

    private double[] lSolve(double D[][], double y[]) { // solve: D'Dx = D'y

        RealMatrix mD = new Array2DRowRealMatrix(D);
        RealMatrix mA = mD.transpose().multiply(mD);

        double b[] = mD.transpose().operate(y);

        // use QR? https://commons.apache.org/proper/commons-math/userguide/linear.html
        DecompositionSolver solver = new QRDecomposition(mA).getSolver();

        RealVector c = new ArrayRealVector(b, false);
        return solver.solve(c).toArray();
    }

    private double w[];
    private double A[];
    private double phi[];
    private double C;

    public void fit4n(double x[], double y[], double w0[]) {

        int nx = x.length, nw = w0.length;
        if (y.length != nx) { throw new RuntimeException("nx != ny"); }

        w = new double[nw];
        for (int i = 0; i < nw; ++i) { w[i] = w0[i]; }

        double D0[][] = new double[nx][2 * nw + 1];

        for (int j = 0; j < nw; ++j) {
            for (int i = 0; i < nx; ++i) {

                double t = w[j] * x[i];
                D0[i][2 * j    ] = Math.cos(t);
                D0[i][2 * j + 1] = Math.sin(t);
            }
        }

        for (int i = 0; i < nx; ++i) { D0[i][2 * nw] = 1.; }

        double est0[] = lSolve(D0, y);
        int nEst = est0.length + nw;

        double estPrev[] = new double[nEst];
        for (int i = 0; i < nEst; ++i) {
            estPrev[i] = (i < est0.length ? est0[i] : 0.);
        }

        for (int it = 0; it < MAXITER; ++it) {

            double D[][] = new double[nx][3 * nw + 1];

            for (int j = 0; j < nw; ++j) {

                double a = estPrev[2 * j    ];
                double b = estPrev[2 * j + 1];

                for (int i = 0; i < nx; ++i) {

                    double t = w[j] * x[i];
                    D[i][2 * j    ] = Math.cos(t);
                    D[i][2 * j + 1] = Math.sin(t);
                    D[i][2 * nw + 1 + j] = 
                        x[i] * ( -a * Math.sin(t) + b * Math.cos(t) );
                }
            }

            for (int i = 0; i < nx; ++i) { D[i][2 * nw] = 1.; }

            double est[] = lSolve(D, y);

            double eps = 0.;
            for (int i = 0; i < est.length; ++i) {
                double tmp = est[i] - estPrev[i];
                eps += tmp * tmp;
            }

            for (int i = 0; i < nw; ++i) {
                w[i] += est[2 * nw + 1 + i];
            }

            for (int i = 0; i < est.length; ++i) { estPrev[i] = est[i]; }

            if (Math.sqrt(eps) < EPS) { break; }
        } // it

        A   = new double[nw];
        phi = new double[nw];

        for (int j = 0; j < nw; ++j) {

            double a = estPrev[2 * j    ];
            double b = estPrev[2 * j + 1];
            double tmp[] = toPolar(a, -b);
            A[j]   = tmp[0];
            phi[j] = tmp[1];
        }

        C = estPrev[2 * nw];
    }

    public void fit4n_wMult(double x[], double y[], double w0, int nw) {

        int nx = x.length;
        if (y.length != nx) { throw new RuntimeException("nx != ny"); }

        double ww = w0;

        double D0[][] = new double[nx][2 * nw + 1];

        for (int j = 0; j < nw; ++j) {
            for (int i = 0; i < nx; ++i) {

                double t = ww * (j + 1) * x[i];
                D0[i][2 * j    ] = Math.cos(t);
                D0[i][2 * j + 1] = Math.sin(t);
            }
        }

        for (int i = 0; i < nx; ++i) { D0[i][2 * nw] = 1.; }

        double est0[] = lSolve(D0, y);
        int nEst = est0.length + 1;

        double estPrev[] = new double[nEst];
        for (int i = 0; i < nEst - 1; ++i) { estPrev[i] = est0[i]; }
        estPrev[nEst - 1] = 0.; // just in case..

        for (int it = 0; it < MAXITER; ++it) {

            double D[][] = new double[nx][2 * nw + 2];

            for (int j = 0; j < nw; ++j) {
                for (int i = 0; i < nx; ++i) {

                    double t = ww * (j + 1) * x[i];
                    D[i][2 * j    ] = Math.cos(t);
                    D[i][2 * j + 1] = Math.sin(t);

                    if (j == 0) {

                        double a = estPrev[0], b = estPrev[1];
                        D[i][2 * nw + 1] = 
                            x[i] * ( -a * Math.sin(t) + b * Math.cos(t) );
                    }
                }
            }

            for (int i = 0; i < nx; ++i) { D[i][2 * nw] = 1.; }

            double est[] = lSolve(D, y);

            double eps = 0.;
            for (int i = 0; i < est.length; ++i) {
                double tmp = est[i] - estPrev[i];
                eps += tmp * tmp;
            }

            ww += est[2 * nw + 1];

            for (int i = 0; i < est.length; ++i) { estPrev[i] = est[i]; }

            if (Math.sqrt(eps) < EPS) { break; }
        } // it

        w = new double[nw];
        for (int i = 0; i < nw; ++i) { w[i] = ww * (i + 1); }

        A   = new double[nw];
        phi = new double[nw];

        for (int j = 0; j < nw; ++j) {

            double a = estPrev[2 * j    ];
            double b = estPrev[2 * j + 1];
            double tmp[] = toPolar(a, -b);
            A[j]   = tmp[0];
            phi[j] = tmp[1];
        }

        C = estPrev[2 * nw];
    }

    private static double cmp(double x[], double y[]) {

        int nx = x.length;
        if (y.length != nx) { throw new RuntimeException("nx != ny"); }

        double d = 0.;
        for (int i = 0; i < nx; ++i) {
            d = Math.max(d, Math.abs(x[i] - y[i]));
        }
        return d;
    }

    private static double[] quantization(double x[], double h) {

        int nx = x.length;
        double v[] = new double[nx];

        double c = 1. / h;
        for (int i = 0; i < nx; ++i) { v[i] = h * Math.round(c * x[i]); }

        return v;
    }


    private static void test_1(double h_q) {

        double h = 1. / F_Q;

        double f_0 = 50.; // [Hz]
        double w_0 = 2. * Math.PI * f_0;

        int nT = 10;
        double T = nT / f_0;

        int N = (int) (T * F_Q);
        double x[] = new double[N];
        for (int i = 0; i < N; ++i) { x[i] = h * i; }

        int nw = 35;

        double wSig[] = Arrays.copyOf(Example.w, nw);
        // some deviations..
        wSig[ 5] += 3.;
        wSig[10] -= 5.;

        double ASig[]   = Arrays.copyOf(Example.A,   nw);
        double phiSig[] = Arrays.copyOf(Example.phi, nw);

        double s[] = Signal.generate(x, ASig, wSig, phiSig, Example.C);

        double y[];
        if (h_q > 1.e-7) {
            y = quantization(s, h_q);
        } else {
            y = s;
        }

        double w0[] = new double[nw];
        for (int i = 0; i < nw; ++i) { w0[i] = w_0 * (i + 1); }

        SineFit sf = new SineFit();
        sf.fit4n(x, y, w0);

        System.out.println("");
        System.out.println(">> test 1, nw = " + nw + ", h_q = " + h_q);
        System.out.println("w:   " + cmp(sf.w,   wSig));
        System.out.println("A:   " + cmp(sf.A,   ASig));
        System.out.println("phi: " + cmp(sf.phi, phiSig));
        System.out.println("C:   " + Math.abs(sf.C));
    }

    private static void test_2(double h_q) {

        double h = 1. / F_Q;

        double f_0 = 50.; // [Hz]
        double w_0 = 2. * Math.PI * f_0;

        int nT = 10;
        double T = nT / f_0;

        int N = (int) (T * F_Q);
        double x[] = new double[N];
        for (int i = 0; i < N; ++i) { x[i] = h * i; }

        int nw = 50;

        double wSig[]   = Arrays.copyOf(Example.w,   nw);
        double ASig[]   = Arrays.copyOf(Example.A,   nw);
        double phiSig[] = Arrays.copyOf(Example.phi, nw);

        double s[] = Signal.generate(x, ASig, wSig, phiSig, Example.C);

        double y[];
        if (h_q > 1.e-7) {
            y = quantization(s, h_q);
        } else {
            y = s;
        }

        SineFit sf = new SineFit();
        sf.fit4n_wMult(x, y, w_0, nw);

        System.out.println("");
        System.out.println(">> test 2, nw = " + nw + ", h_q = " + h_q);
        System.out.println("w:   " + cmp(sf.w,   wSig));
        System.out.println("A:   " + cmp(sf.A,   ASig));
        System.out.println("phi: " + cmp(sf.phi, phiSig));
        System.out.println("C:   " + Math.abs(sf.C));
    }

    public static void main(String args[]) {

        test_1(0.);
        test_1(0.001);

        test_2(0.);
        test_2(0.001);
    }
}
