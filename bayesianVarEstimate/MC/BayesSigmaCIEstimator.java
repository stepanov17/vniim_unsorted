
import java.util.Random;


public class BayesSigmaCIEstimator {

    private static double[] getSigmaRange(double step, double maxSigma) {

        int n = (int) (maxSigma / step) + 1;
        double r[] = new double[n];
        for (int i = 0; i < n; ++i) {
            r[i] = i * step;
        }

        r[0] = 1.e-10; // avoid div. by zero

        return r;
    }

    private static double gamma(double x) {

        double p[] = {
             1.000000000190015, 76.18009172947146, -86.50532032941677, 24.01409824083091,
            -1.231739572450155, 1.208650973866179e-3, -5.395239384953e-6};

        double res = p[0];
        for (int i = 1; i < p.length; ++i) {
            res += p[i] / (x + (double) i);
        }

        res *= (Math.sqrt(2. * Math.PI) / x);
        res *= Math.pow(x + 5.5, x + 0.5);
        res *= Math.exp(-(x + 5.5));

        return res;
    }

    private final Random rnd = new Random(100500L);

    private final double p0 = 0.95;

    private final double maxSigma;
    private final double dSigma;

    private final double sigma[];
    private final double pdf1[], pdf2[];
    //private double E1 = Double.NaN, E2 = Double.NaN;

    private final int nSigma;

    private final int n1, n2;
    private final double sigma1, sigma2;

//    double minLeft = 1.e10, maxRight = 0.;
    double maxRight = 0.;

    public BayesSigmaCIEstimator(int    n1,
                                 double sigma1,
                                 int    n2,
                                 double sigma2,
                                 double maxSigma,
                                 double dSigma) {
        this.n1     = n1;
        this.sigma1 = sigma1;
        this.n2     = n2;
        this.sigma2 = sigma2;

        this.maxSigma = maxSigma;
        this.dSigma   = dSigma;

        sigma = getSigmaRange(dSigma, maxSigma);
        nSigma = sigma.length;

        pdf1 = new double[nSigma];
        pdf2 = new double[nSigma];
    }

    private double[] getCenteredSample(int n, double s) {

        double res[] = new double[n];
        double mean = 0.;
        for (int i = 0; i < n; ++i) {
            double v = s * rnd.nextGaussian();
            res[i] = v;
            mean += v;
        }

        mean /= n;

        for (int i = 0; i < n; ++i) { res[i] -= mean; }

        return res;
    }

    private void recalculatePDFs(double sample1[], double sample2[]) {

        if ((sample1.length != n1) || (sample2.length != n2)) {
            throw new RuntimeException("invalid sample length");
        }

        double S = 0.;

        for (double v: sample1) { S += (v * v); }
        double S1 = 0.5 * S;

        for (double v: sample2) { S += (v * v); }
        double S2 = 0.5 * S;

        double C1 = (2. * Math.pow(S1, 0.5 * (n1 - 1))) / gamma(0.5 * (n1 - 1));
        double C2 = (2. * Math.pow(S2, n2 - 1)) / gamma(n2 - 1);

        //E1 = Math.sqrt(S1) * gamma(0.5 * n1 - 1.) / gamma(0.5 * (n1 - 1));
        //E2 = Math.sqrt(S2) * gamma(n2 - 1.5) / gamma(n2 - 1);

        int p1 = n1;
        int p2 = 2 * (n2 - 1) + 1;

        for (int i = 0; i < nSigma; ++i) {

            double s = sigma[i];
            pdf1[i] = C1 * Math.exp(-S1 / (s * s)) * Math.pow(s, -p1);
            pdf2[i] = C2 * Math.exp(-S2 / (s * s)) * Math.pow(s, -p2);
        }
    }

    private double coverageIntervalLength(double pdf[], int i0) {

        double I = 0.;

        int i1 = -1;
        double f1 = pdf[i0], f2;
        for (int i = i0 + 1; i < nSigma; ++i) {

            f2 = pdf[i];
            I += 0.5 * (f1 + f2) * dSigma;
            if (I >= p0) {
                i1 = i;
                // debug
                //System.out.println(i1 * dSigma);
                break;
            }
            f1 = f2;
        }

        if (I < p0) {
            return -1.;
        }

        //System.out.println("I ~ " + I); // debug
        if (Math.abs(I - p0) > 1.e-3) { // should not exceed 0.1%
            System.err.println("too rough estimate, I ~ " + I);
        }

        return (i1 - i0) * dSigma;
    }

    private int coverageIntervalStartIndex(double pdf[], double eps) {

        double I = 0.;

        double f1 = pdf[0];
        for (int i = 1; i < nSigma; ++i) {

            double f2 = pdf[i];
            I += 0.5 * (f1 + f2) * dSigma;
            if (I >= eps) {
                return i;
            }
            f1 = f2;
        }

        return 0;
    }

    // return a pair [start_point, length]
    private double[] coverageInterval(double pdf[]) {

        double ci[] = new double[]{-1., -1.};

        int i0 = coverageIntervalStartIndex(pdf, 1.e-5);
        //System.out.println("start from " + i0);

        double L0 = 10. * maxSigma;

        int nFound = 0;
        for (int i = i0; i < nSigma; ++i) {
            double L = coverageIntervalLength(pdf, i);
            if (L < 0) {
                break;
            }
            if (L < L0) {
                L0 = L;
                i0 = i;
                ++nFound;
            }
        }

        if (nFound > 0) {
            //ci[0] = i0 * dSigma;
            ci[0] = (i0 + 1) * dSigma; //+1??
            ci[1] = L0;
        }

        return ci;
    }

    private int iteration() {

        double x1[] = getCenteredSample(n1, sigma1);
        double x2[] = getCenteredSample(n2, sigma2);

        recalculatePDFs(x1, x2);

        double ci1[] = coverageInterval(pdf1);
        if (ci1[1] < 0.) { return -1; }

        double ci2[] = coverageInterval(pdf2);
        if (ci2[1] < 0.) { return -1; }

        double l1 = ci1[1], l2 = ci2[1];

        maxRight = Math.max(maxRight, ci1[0] + l1);
        maxRight = Math.max(maxRight, ci2[0] + l2);

        if (l1 > l2) { return 1; } // todo: other checks
        return 0;
    }

    private void MC(int nSim) {

        double P = 0.;
        int n = 0;
        for (int sim = 1; sim <= nSim; ++sim) {

            int v = iteration();
            if (v > -1) {
                P += v;
                ++n;
            }

            System.out.printf(">> %d >> %.3f\n", sim, P / n);
        }

        double u = 1. - (n + 0.) / nSim;
        System.out.printf("undefined: %d (%.2f%%)\n", (nSim - n), 100. * u);
    }

    private void printMaxRight() {
        System.out.printf(">> maxRight = %.5f\n", maxRight);
    }

    public static void main(String[] args) {

        int n = 10;
        BayesSigmaCIEstimator e = new BayesSigmaCIEstimator(
            n, 1., n, 1., 20., 2.e-4);

        e.MC(100_000);

        System.out.println("");

        e.printMaxRight();
    }
}
