
import java.util.Arrays;
import java.util.Random;

/**
 * @author astepanov
 */
public class VarEqualityTests {

    public static enum TEST{VAR_RATIO, RANGE_RATIO};

    private final TEST testType;

    public VarEqualityTests(TEST testType) {  this.testType = testType; }

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

    private double getVarRatio(double x1[], double x2[]) {

        int n1 = x1.length, n2 = x2.length;

        double m1 = 0., m2 = 0.;

        for (double v: x1) { m1 += v; }
        m1 /= n1;

        for (double v: x2) { m2 += v; }
        m2 /= n2;

        double s1 = 0., s2 = 0.;

        for (double v: x1) { s1 += (v - m1) * (v - m1); }
        for (double v: x2) { s2 += (v - m2) * (v - m2); }

        s1 /= (n1 - 1.);
        s2 /= (n2 - 1.);

        double r = s1 / s2;
        return (s1 > s2) ? r : 1. / r;
    }

    private double getRangeRatio(double x1[], double x2[]) {

        int n1 = x1.length, n2 = x2.length;

        Arrays.sort(x1);
        Arrays.sort(x2);

        double r1 = x1[n1 - 1] - x1[0];
        double r2 = x2[n2 - 1] - x2[0];

        double r = r1 / r2;
        return (r1 > r2) ? r : 1. / r;
    }

    private double getR(double x1[], double x2[], TEST type) {

        if (type == TEST.VAR_RATIO) {
            return getVarRatio(x1, x2);
        }

        return getRangeRatio(x1, x2);
    }

    private double[] getCriticalVals(int n1, int n2, double p, int nSim) {

        double r[] = new double[nSim];
        for (int i = 0; i < nSim; ++i) {

            double x1[] = getTSPSample(n1, p);
            double x2[] = getTSPSample(n2, p);

            r[i] = getR(x1, x2, testType);
        }

        Arrays.sort(r);

        int i1 = (int) (0.025 * nSim), i2 = (int) (0.975 * nSim);

        double res[] = new double[2];
        res[0] = r[i1 - 1];
        res[1] = r[i2 - 1];
        return res;
    }

    private double getP(int n1, int n2, double p, double c, double C1, double C2, int nSim) {

        double P = nSim;

        for (int i = 0; i < nSim; ++i) {

            double x1[] = getTSPSample(n1, p);
            double x2[] = getTSPSample(n2, p);

            if (Math.abs(c - 1.) > 1.e-8) {
                for (int j = 0; j < n2; ++j) { x2[j] *= c; }
            }

            double r = getR(x1, x2, testType);
            if ((r < C1) || (r > C2)) { --P; }
        }

        return P / nSim;
    }

    private static double roundTo(double x, double h) {

        double tmp = 1. / h;
        return h * Math.round(tmp * x);
    }

    public static void main(String[] args) {

        VarEqualityTests calculator =
                new VarEqualityTests(VarEqualityTests.TEST.VAR_RATIO);

        int N[][] = {
            {10,  5},
            {10, 10},
            {15,  5},
            {15, 10},
            {15, 15},
            {20,  5},
            {20, 10},
            {20, 15},
            {20, 20}
        };

        int nSim = 1_000_000, nAvg = 200;

        double p = 1.;

        for (int n[]: N) {

            int n1 = n[0], n2 = n[1];

            double c1 = 0., c2 = 0.;
            for (int i = 0; i < nAvg; ++i) {
                double cv[] = calculator.getCriticalVals(n1, n2, p, nSim);
                c1 += cv[0];
                c2 += cv[1];
            }
            c1 /= nAvg;
            c2 /= nAvg;

            c1 = roundTo(c1, 1.e-3);
            c2 = roundTo(c2, 1.e-2);

            // check (must be ~0.95)
            double P = calculator.getP(n1, n2, p, 1., c1, c2, 5 * nSim);

            System.out.printf("%d\t%d\t%.3f\t%.2f\t%.3f\n", n1, n2, c1, c2, P);
        }
    }
}
