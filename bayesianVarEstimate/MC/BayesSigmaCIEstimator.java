import java.util.Arrays;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadLocalRandom;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;


public class BayesSigmaCIEstimator {

    private final static int MC_TIMEOUT_SEC = 3 * 3600;

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

    private final double p0 = 0.95;

    private final double maxSigma;
    private final double dSigma;
    private final double sigma[];
    private final int nSigma;

    private final int n1, n2;
    private final double sigma1, sigma2;


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
    }

    private double[] getCenteredSample(int n, double s) {

        ThreadLocalRandom rnd = ThreadLocalRandom.current();

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

    private double iteration() {

        double sample1[] = getCenteredSample(n1, sigma1);
        double sample2[] = getCenteredSample(n2, sigma2);

        double S = 0.;

        for (double v: sample1) { S += (v * v); }
        double S1 = 0.5 * S;

        for (double v: sample2) { S += (v * v); }
        double S2 = 0.5 * S;

        double C1 = (2. * Math.pow(S1, 0.5 * (n1 - 1))) / gamma(0.5 * (n1 - 1));
        double C2 = (2. * Math.pow(S2, n2 - 1)) / gamma(n2 - 1);

        //double E1 = Math.sqrt(S1) * gamma(0.5 * n1 - 1.) / gamma(0.5 * (n1 - 1));
        //double E2 = Math.sqrt(S2) * gamma(n2 - 1.5) / gamma(n2 - 1);

        int p1 = n1;
        int p2 = 2 * (n2 - 1) + 1;

        double pdf1[] = new double[nSigma];
        double pdf2[] = new double[nSigma];

        for (int i = 0; i < nSigma; ++i) {

            double s = sigma[i];
            pdf1[i] = C1 * Math.exp(-S1 / (s * s)) * Math.pow(s, -p1);
            pdf2[i] = C2 * Math.exp(-S2 / (s * s)) * Math.pow(s, -p2);
        }

        double ci1[] = coverageInterval(pdf1);
        double ci2[] = coverageInterval(pdf2);

        double l1 = ci1[1], l2 = ci2[1];
        return l2 / l1; // todo: other checks
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

    private double MC(int nSim, int nThreads, boolean calculateCritValue, double C) throws InterruptedException {

        ThreadPoolExecutor executor =
            (ThreadPoolExecutor) Executors.newFixedThreadPool(nThreads);

        double r[] = new double[nSim];

        Arrays.fill(r, -1.);

        for (int sim = 0; sim < nSim; ++sim) {

            int i = sim;

            executor.execute(

                () -> {

                    r[i] = iteration();

                    if ((i > 0) && ((i + 1) % 1000 == 0)) {
                        System.out.println(">> " + (i + 1));
                    }
                });
        }

        executor.shutdown();

        Thread.sleep(100);
        if (!executor.awaitTermination(MC_TIMEOUT_SEC, TimeUnit.SECONDS)) {
            throw new RuntimeException("the calculations were not completed");
        }

        int undefined = 0;
        double res;

        if (calculateCritValue) {

            Arrays.parallelSort(r);

            for (double v: r) {
                if (v < 0.) { ++undefined; }
                else { break; }
            }

            double distr[] = Arrays.copyOfRange(r, undefined, nSim);
            int i0 = (int) (p0 * distr.length);
            res = distr[i0 - 1];

        } else {

            int s = 0;
            for (double v: r) {
                if (v < 0.) { ++undefined; }
                else if (v <= C) { ++s; }
            }
            res = s;
            res /= (nSim - undefined);
        }

        if (undefined > 0) { System.err.println("undefined: " + undefined); }

        return res;
    }

    public double calculateCriticalValue(int nSim, int nThreads) throws InterruptedException {

        return MC(nSim, nThreads, true, -1.);
    }

    public double calculateP(double cv, int nSim, int nThreads) throws InterruptedException {

        return MC(nSim, nThreads, false, cv);
    }


    public static void main(String args[]) throws InterruptedException {

        int n = 10, nSim = 100_000;

        int nThreads = Math.max(1, Runtime.getRuntime().availableProcessors() / 2);
        System.out.println("nThreads = " + nThreads);

        BayesSigmaCIEstimator e = new BayesSigmaCIEstimator(
            n, 1., n, 1., 15., 2.e-4);

        double cv = e.calculateCriticalValue(nSim, nThreads);
        System.out.printf("%.3f\n\n", cv);

        // check: ~ p0
        double p = e.calculateP(cv, nSim / 2, nThreads);
        System.out.printf("%.3f\n", p);
    }
}
