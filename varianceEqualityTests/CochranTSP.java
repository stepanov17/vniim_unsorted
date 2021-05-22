import java.util.Arrays;
import java.util.Random;


public class CochranTSP {

    private final static Random RND = new Random(/*123456L*/);

    private final static double P0 = 0.95;

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

    double getS2(double x[]) {

        int nx = x.length;

        double m = 0., s = 0.;
        for (double v: x) { m += v; }
        m /= nx;

        for (double v: x) { s += (v - m) * (v - m); }

        return s / (nx - 1.);
    }

    double getCochran(double s2[]) {

        if (s2.length < 2) {
            throw new IllegalArgumentException("wrong number of samples");
        }

        double sum = s2[0];
        double maxs = s2[0];
        for (int i = 1; i < s2.length; ++i) {
            double v = s2[i];
            if (v > maxs) { maxs = v; }
            sum += v;
        }
        return maxs / sum;
    }

    private double getCriticalVal(int nx, int nSamples, double p, int nSim) {

        double r[] = new double[nSim];

        for (int i = 0; i < nSim; ++i) {

            double s2[] = new double[nSamples];
            for (int j = 0; j < nSamples; ++j) {
                s2[j] = getS2(getTSPSample(nx, p));
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
                s2[j] = getS2(x);
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

        CochranTSP calculator = new CochranTSP();

        int nSim = 1_000_000, nAvg = 20;
        double p = 2;

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
