
import java.util.Arrays;
import java.util.Random;


public class CochranProximity {

    private final double dR;

    public static enum TEST {
        VAR_RATIO,    // Cochran
        RANGE_RATIO,  // Bliss, Cochran, Tukey
    };

    private final static Random RND = new Random();

    private final static double P0 = 0.95;

    private final TEST testType;

    public CochranProximity(TEST testType, double dR) {
        this.testType = testType;
        this.dR = dR;
    }

    // p > 0 => TSP
    // p < 0 => normal
    private double[] getTSPSample(int n, double p) {

        double e[] = new double[n];

        if (p < 0.) {
            for (int i = 0; i < n; ++i) { e[i] = RND.nextGaussian(); }
        } else {

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
        }

        double s = RND.nextDouble() * 2. * dR + 1. - dR;
        for (int j = 0; j < n; ++j) { e[j] *= s; }

        return e;
    }

    private double getS2(double x[]) {

        int nx = x.length;
        if (nx < 3) { throw new IllegalArgumentException("invalid nx: " + nx); }

        double m = 0., s = 0.;
        for (double v: x) { m += v; }
        m /= nx;

        for (double v: x) { s += (v - m) * (v - m); }

        return s / (nx - 1.);
    }

    private double getR(double x[]) {

        int nx = x.length;
        if (nx < 3) { throw new IllegalArgumentException("invalid nx: " + nx); }

        double minx = x[0], maxx = x[0];
        for (int i = 1; i < nx; ++i) {
            double t = x[i];
            if (t < minx) { minx = t; }
            else if (t > maxx) { maxx = t; }
        }

        return maxx - minx;
    }

    private double getS2orR(double x[]) {

        if (testType == TEST.VAR_RATIO) {
            return getS2(x);
        } else {
            return getR(x);
        }
    }

    private double getCochran(double v[]) {

        if (v.length < 2) {
            throw new IllegalArgumentException("wrong number of samples");
        }

        double sum = v[0];
        double maxs = v[0];
        for (int i = 1; i < v.length; ++i) {
            double t = v[i];
            if (t > maxs) { maxs = t; }
            sum += t;
        }
        return maxs / sum;
    }

    private double getCriticalVal(int nx, int nSamples, double p, int nSim) {

        double r[] = new double[nSim];

        for (int i = 0; i < nSim; ++i) {

            double s2[] = new double[nSamples];
            for (int j = 0; j < nSamples; ++j) {
                s2[j] = getS2orR(getTSPSample(nx, p));
            }
            r[i] = getCochran(s2);
        }

        Arrays.sort(r);
        int i0 = (int) (P0 * nSim);
        return r[i0 - 1];
    }

    private double getCriticalValAvg(int nx, int nSamples, double p, int nSim, int nAvg) {

        double c = 0.;
        for (int i = 0; i < nAvg; ++i) {
            c += getCriticalVal(nx, nSamples, p, nSim);
        }
        return c / nAvg;
    }

    private double getP(int nx, int nSamples, double p, double c, double critValue, int nSim) {

        double P = nSim;

        for (int i = 0; i < nSim; ++i) {

            double s2[] = new double[nSamples];
            for (int j = 0; j < nSamples; ++j) {

                double x[] = getTSPSample(nx, p);

                if ((j == 0) && (Math.abs(c - 1.) > 1.e-8)) {
                    for (int k = 0; k < nx; ++k) { x[0] *= c; }
                }
                s2[j] = getS2orR(x);
            }
            double r = getCochran(s2);
            if (r > critValue) { --P; }
        }

        return P / nSim;
    }

    private double getPAvg(int nx, int nSamples, double p, double c, double critValue, int nSim, int nAvg) {

        double P = 0.;
        for (int i = 0; i < nAvg; ++i) {
            P += getP(nx, nSamples, p, c, critValue, nSim);
        }
        return P / nAvg;
    }


    public static void main(String args[]) {

        double dR = 0.10;

        CochranProximity calc = new CochranProximity(TEST.VAR_RATIO, dR);

        int nSamples = 5, nx = 10;
        double c = calc.getCriticalValAvg(nx, nSamples, -1., 1_000_000, 100);
        System.out.printf("%.3f\n", c);

        // check: should be 0.95
        double P = calc.getPAvg(nx, nSamples, -1., 1., c, 1_000_000, 10);
        System.out.printf("%.3f\n", P);
    }
}
