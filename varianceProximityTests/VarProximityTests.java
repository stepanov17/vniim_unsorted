
import java.util.Arrays;
import java.util.Random;


/**
 * @author astepanov
 */
public class VarProximityTests {

    private final double dR;

    //  I) H0 : sigma_2 / sigma_1 ~ U(1 - dR, 1 + dR)
    // II) H0 : sigma_1, sigma_2 = sigma_0 * d, d ~ U(1 - dR, 1 + dR)
    public static enum HYPOTHESIS{I, II};

    public static enum TEST{VAR_RATIO, RANGE_RATIO};

    private final HYPOTHESIS hypothesis;
    private final TEST testType;

    public VarProximityTests(double dR, HYPOTHESIS hypothesis, TEST testType) {

        this.dR         = dR;
        this.hypothesis = hypothesis;
        this.testType   = testType;
    }

    private final static Random RND = new Random(/*123456L*/);

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

    private double getRange(double x[]) {

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

    private double getR(double x1[], double x2[]) {

        if (testType == TEST.RANGE_RATIO) {
            return getRange(x1) / getRange(x2);
        } else {
            return getS2(x1) / getS2(x2);
        }
    }


    private double[] getCriticalVals(int n1, int n2, double p, int nSim) {

        double r[] = new double[nSim];
        for (int i = 0; i < nSim; ++i) {

            double x1[] = getTSPSample(n1, p);
            double x2[] = getTSPSample(n2, p);

            if (hypothesis == HYPOTHESIS.II) {

                double s1 = RND.nextDouble() * 2. * dR + 1. - dR;
                for (int j = 0; j < n1; ++j) { x1[j] *= s1; }
            }

            // I, II
            double s2 = RND.nextDouble() * 2. * dR + 1. - dR;
            for (int j = 0; j < n2; ++j) { x2[j] *= s2; }

            r[i] = getR(x1, x2);
        }

        Arrays.sort(r);

        int i1 = (int) (0.025 * nSim), i2 = (int) (0.975 * nSim);

        double res[] = new double[2];
        res[0] = r[i1 - 1];
        res[1] = r[i2 - 1];
        return res;
    }

    private double[] getCriticalValsAvg(int n1, int n2, double p, int nSim, int nAvg) {

        double c1 = 0., c2 = 0.;
        for (int i = 0; i < nAvg; ++i) {
            double c[] = getCriticalVals(n1, n2, p, nSim);
            c1 += c[0];
            c2 += c[1];
        }
        return new double[]{c1 / nAvg, c2 / nAvg};
    }

    private double getP(int n1, int n2, double p, double c, double C1, double C2, int nSim) {

        double P = nSim;

        for (int i = 0; i < nSim; ++i) {

            double x1[] = getTSPSample(n1, p);
            double x2[] = getTSPSample(n2, p);

            if (hypothesis == HYPOTHESIS.II) {

                double s1 = RND.nextDouble() * 2. * dR + 1. - dR;
                for (int j = 0; j < n1; ++j) { x1[j] *= s1; }
            }

            // I, II
            double s2 = c * (RND.nextDouble() * 2. * dR + 1. - dR);
            for (int j = 0; j < n2; ++j) { x2[j] *= s2; }

            double r = getR(x1, x2);
            if ((r < C1) || (r > C2)) { --P; }
        }

        return P / nSim;
    }

    private double getPAvg(int n1, int n2, double p, double c, double C1, double C2, int nSim, int nAvg) {

        double P = 0.;
        for (int i = 0; i < nAvg; ++i) {
            P += getP(n1, n2, p, c, C1, C2, nSim);
        }
        return P / nAvg;
    }


    public static void main(String args[]) {

        double dR = 0.1; // sigma +/-10%

        VarProximityTests calculator = new VarProximityTests(
                dR,
                HYPOTHESIS.II,
                TEST.VAR_RATIO);

        int N[][] = {
            {5, 5},
            {10, 10},
            {20, 20}
        };

        int nSim = 1_000_000, nAvg = 500;

        double p = 1.;

        for (int n[]: N) {

            int n1 = n[0], n2 = n[1];

            double c[] = calculator.getCriticalValsAvg(n1, n2, p, nSim, nAvg);
            double c1 = c[0], c2 = c[1];

            // check (must be ~0.95)
            double P = calculator.getPAvg(n1, n2, p, 1., c1, c2, nSim, 10);

            System.out.printf("%d\t%d\t%.3f\t%.3f\t%.3f\n", n1, n2, c1, c2, P);
        }
    }
}
