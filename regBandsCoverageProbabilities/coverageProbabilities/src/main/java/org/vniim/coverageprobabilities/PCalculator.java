
package org.vniim.coverageprobabilities;

import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;


public class PCalculator {

    private final int n;
    private final int n0;
    private final double a;

    private final double rho;
    private final double acos_rho;


    public PCalculator(int n, int n0, double a) {

        this.n  = n;
        this.n0 = n0;
        this.a  = a;

        rho = calculateRho();
        acos_rho = Math.acos(rho);
    }

    private double calculateRho() {

        double r[] = new double[n0];
        r[0] = 1.;
        for (int i = 1; i < n0; ++i) { r[i] = a * r[i - 1]; }

        double vals[][] = new double[n0][n0];
        for (int i = 0; i < n0; ++i) {
            int k = 0;
            for (int j = i; j < n0; ++j) {
                vals[i][j] = r[k];
                vals[j][i] = r[k];
                ++k;
            }
        }

        RealMatrix V = MatrixUtils.createRealMatrix(vals);
        RealMatrix VInv = new LUDecomposition(V).getSolver().getInverse();

        vals = new double[n0][2];
        for (int i = 0; i < n0; ++i) {
            vals[i][0] = 1.;
            vals[i][1] = i;
        }

        RealMatrix A = MatrixUtils.createRealMatrix(vals);

        RealMatrix Prod = A.transpose().multiply(VInv).multiply(A);
        RealMatrix T = new LUDecomposition(Prod).getSolver().getInverse();

        vals = new double[][]{{1., 0.}};
        RealMatrix v1 = MatrixUtils.createRealMatrix(vals);

        vals = new double[][]{{1., n - 1.}};
        RealMatrix v2 = MatrixUtils.createRealMatrix(vals);

        double num = -v1.multiply(T).multiply(v2.transpose()).getEntry(0, 0);
        double den = Math.sqrt(
            v1.multiply(T).multiply(v1.transpose()).getEntry(0, 0) *
            v2.multiply(T).multiply(v2.transpose()).getEntry(0, 0)
        );

        return num / den;
    }

    private double G(double t) {

        double g = 1.;
        if ((t < 1.) && (t >= Math.sqrt(0.5 + 0.5 * rho))) {
            g = (acos_rho - 2. * Math.acos(t)) / Math.PI;
        }
        return g;
    }

    private double f2(double x, double r) {

        return G(x / r) * r * Math.exp(-0.5 * r * r);
    }

    private double F(double x) {

        double I1 = 1. - Math.exp(-0.5 * x * x);

        double I2 = 0;
        double r_max = x * Math.sqrt(2. / (1. + rho));
        double h = 1.e-5 * (r_max - x);

        double r = x;
        double f2_prev = f2(x, r);

        while (r <= r_max + 0.1 * h) {

            r += h;
            double f2 = f2(x, r);
            I2 += 0.5 * (f2_prev + f2);
            f2_prev = f2;
        }

        I2 *= h;

        return I1 + I2;
    }

    private static void printP(int n, int n0, double K) {

        System.out.printf("n0 = %4d ", n0);

        for (double a = 0.; a < 0.9501; a += 0.05) {

            PCalculator calc = new PCalculator(n, n0, a);
            System.out.printf(" | %.4f", calc.F(K));
        }
        System.out.println("");
    }

    public static void printTableHead() {

        for (int i = 0; i < 191; ++i) { System.out.print("_"); }
        System.out.println("");
        System.out.print("     alpha");
        for (double a = 0.; a < 0.9501; a += 0.05) {
            System.out.printf(" |  %.2f ", a);
        }
        System.out.println("");
        for (int i = 0; i < 191; ++i) { System.out.print("_"); }
        System.out.println("");
    }


    public static void main(String args[]) {

        double KK[] = {2.1403, 2.4414, 3.0273};
        int PP[] = {90, 95, 99};

        for (int i = 0; i < KK.length; ++i) {

            double K = KK[i];

            System.out.println("");
            System.out.println("K = " + K + " (for coverage probability level of ~" + PP[i] + "%)");
            System.out.println("");
            System.out.println("");

            int nn[] = {50, 100, 200, 400, 600, 800, 1000, 1500, 2000, 3000};

            for (int n: nn) {

                System.out.println("n = " + n);
                System.out.println();

                printTableHead();

                for (int n0 = (int)(0.2 * n); n0 <= (int)(0.8 * n); n0 += (int)(0.1 * n)) {

                    printP(n, n0, K);
                }
                System.out.println();
                System.out.println();
                System.out.println();
            }

            if (i < KK.length - 1) {
                for (int j = 0; j < 191; ++j) { System.out.print("="); }
            }
            System.out.println("");
        }
    }
}
