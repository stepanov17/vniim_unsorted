
import java.util.Arrays;
import java.util.Random;


public class CochranTSP {

    public static enum TEST{
        VAR_RATIO,    // Cochran
        RANGE_RATIO   // Bliss, Cochran, Tukey
    };

    private final static Random RND = new Random();

    private final static double P0 = 0.95;

    private final TEST testType;

    public CochranTSP(TEST testType) {  this.testType = testType; }

    // p > 0 => TSP
    // p < 0 => normal
    private double[] getTSPSample(int n, double p) {

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

    private void calculateTestPow(int nSamples, int nx, double p) {

        double kMax = 5., dk = 5.e-3;

        System.out.println("// " + testType +  ", nSamples = " + nSamples + ", nx = " + nx + ", p = " + p);
        System.out.println("");

        int nSim = 1_000_000, nAvg = 10;

        double c = getCriticalValAvg(nx, nSamples, p, nSim, 2 * nAvg);

        double k = 1., pow = 0.;

        System.out.println("data = [");

        boolean first = true;

        while ((k <= kMax + 0.1 * dk) && (pow <= 0.9995)) {

            pow = 1. - getPAvg(nx, nSamples, p, k, c, nSim, nAvg);
            if (first) {
                // check
                if (Math.abs(1. - P0 - pow) > 2.e-3) {
                    System.err.println("pow != 1 - P0: " + pow);
                }
                first = false;
            }

            System.out.printf("%.3f\t%.3f\n", k, pow);
            k += dk;
        }

        System.out.println("];");
    }


    public static void main(String args[]) {

        CochranTSP calculator = new CochranTSP(TEST.VAR_RATIO);

        int nSim = 1_000_000, nAvg = 20;
        double p = -1.;

        int nSamples = 3, NX[] = {5, 10, 15, 20};

        System.out.println("p = " + p + ", nSamples = " + nSamples);

        for (int nx: NX) {
            double c = calculator.getCriticalValAvg(nx, nSamples, p, nSim, nAvg);
            // check
            double chk = calculator.getPAvg(nx, nSamples, p, 1., c, nSim, 3);
            System.out.printf("nx = %02d:\t%.3f\t(p0 = %.3f)\n", nx, c, chk);
        }
    }
}
