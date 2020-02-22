package org.vniim.arbands;

import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;

public class ForecastCalculator {

    private final double a;
    private final int n, n0;

    private final double R[], U[];

    private final RealMatrix M;
    private final RealMatrix Theta;

    private double[] calcR() { // acf

        double r[] = new double[n];
        r[0] = 1.;
        for (int i = 1; i < n; ++i) { r[i] = a * r[i - 1]; }
        return r;
    }

    public ForecastCalculator(int n, int n0, double a) {

        this.n  = n;
        this.n0 = n0;
        this.a  = a;

        R = calcR();

        double w[][] = new double[n0][n0];
        for (int i = 0; i < n0; ++i) {
            int k = 0;
            for (int j = i; j < n0; ++j) {
                w[i][j] = R[k];
                w[j][i] = R[k];
                ++k;
            }
        }
        RealMatrix W = MatrixUtils.createRealMatrix(w);
        RealMatrix WI = new LUDecomposition(W).getSolver().getInverse();

        double v[][] = new double[n0][2];
        for (int i = 0; i < n0; ++i) {
            v[i][0] = 1.; v[i][1] = i;
        }
        RealMatrix V = MatrixUtils.createRealMatrix(v);
        RealMatrix T = V.transpose().multiply(WI).multiply(V);
        Theta = new LUDecomposition(T).getSolver().getInverse();
        M = Theta.multiply(V.transpose()).multiply(WI);

        double
                theta_11 = Theta.getEntry(0, 0),
                theta_12 = Theta.getEntry(0, 1),
                theta_22 = Theta.getEntry(1, 1);

        U = new double[n];
        for (int i = 0; i < n; ++i) {
            U[i] = Math.sqrt(theta_11 + 2 * i * theta_12 + i * i * theta_22);
        }
    }

    public double[] getU() { return U; }

    public double[] getForecast(double data[]) {

        if (data.length != n) { throw new RuntimeException("invalid data size"); }
        RealMatrix S = MatrixUtils.createRealMatrix(n0, 1);
        for (int i = 0; i < n0; ++i) { S.setEntry(i, 0, data[i]); }
        RealMatrix C = M.multiply(S);
        double a = C.getEntry(0, 0), b = C.getEntry(1, 0);

        double forecast[] = new double[n];
        for (int i = 0; i < n; ++i) { forecast[i] = a + b * i; }
        return forecast;
    }

    public void printTheta() {
        for (int i = 0; i < Theta.getRowDimension(); ++i) {
            for (int j = 0; j < Theta.getColumnDimension(); ++j) {
                System.out.print(Theta.getEntry(i, j) + "  ");
            }
            System.out.println("");
        }
        System.out.println(M.getEntry(0, 0));
    }
}
