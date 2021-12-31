
import java.util.Arrays;

public class BayesSigmaCI {

    private final double samples[][] = {

        {4.6156, 3.4793, 3.7839, 3.0842},
        {2.9247, 3.4047, 3.7965, 3.8789},
        {2.9844, 2.8869, 0.9911, 1.6435},
        {3.6131, 3.2003, 3.0293, 1.3871},
        {3.7788, 3.6189, 5.0602, 3.6211},
        {2.5377, 4.1953, 4.1078, 3.3730},
        {3.5031, 1.5954, 3.6658, 3.6223},
        {3.5038, 1.4432, 2.3405, 1.4001},
        {1.4236, 1.4021, 3.4086, 1.9935},
        {4.1279, 4.0914, 2.9162, 2.5567},
        {2.5955, 3.6656, 2.8175, 1.9649},
        {2.7929, 3.9223, 2.6684, 2.1350},
        {2.0403, 1.8503, 1.1314, 2.4318},
        {1.2347, 1.4374, 1.5545, 1.6394},
        {5.6078, 4.4117, 3.5977, 2.5501},
        {4.3074, 5.7166, 3.7637, 5.0008},
        {1.9076, 2.8437, 2.7861, 0.7798},
        {4.5797, 3.7341, 2.2552, 3.7867},
        {1.4697, 2.0862, 1.7895, 1.6111},
        {3.2589, 3.1021, 2.7143, 4.7598},
        {3.4456, 3.8681, 5.6861, 5.7957},
        {3.4702, 5.6231, 3.2774, 3.0398},
        {1.3291, 3.4120, 3.0197, 4.3298},
        {1.3899, 3.6880, 4.6593, 1.2786},
        {3.0037, 4.7318, 2.5069, 2.2075}
    };

    private final int nSamples = samples.length;

    private double S[];
    private double C[];
    private int P[];

    private double E[];

    private final double p0 = 0.95;


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

    private void centerSamples() {

        for (int i = 0; i < nSamples; ++i) {

            double m = Arrays.stream(samples[i]).average().getAsDouble();
            samples[i] = Arrays.stream(samples[i]).map(v -> v - m).toArray();
        }
    }

    private void calculatePdfParams() {

        centerSamples();

        S = new double[nSamples];
        C = new double[nSamples];
        P = new int[nSamples];

        E = new double[nSamples];

        double s = 0.;

        for (int i = 0; i < nSamples; ++i) {

            s += 0.5 * Arrays.stream(samples[i]).map(d -> d * d).sum();
            S[i] = s;

            double p = 0.5 * (i + 1.) * (samples[i].length - 1);
            C[i] = 2. * Math.pow(S[i], p) / gamma(p);

            E[i] = Math.sqrt(s) * (gamma(p - 0.5) / gamma(p));

            P[i] = (i + 1) * (samples[i].length - 1) + 1;
        }
    }

    private double[] pdf(int i, double dSigma, int n) {

        if ((i < 0) || (i >= nSamples)) { 
            throw new IllegalArgumentException("invalid index"); }

        double pdf[] = new double[n];
        for (int j = 1; j < n; ++j) {

            double x = j * dSigma;
            pdf[j] = C[i] * Math.pow(x, -P[i]) * Math.exp(-S[i] / (x * x));
            if (Double.isNaN(pdf[j])) {
                throw new RuntimeException("NaN in pdf");
            }
        }

        pdf[0] = 0.;

        return pdf;
    }

    private int coverageIntervalStartIndex(double pdf[], double dSigma) {

        double I = 0.;

        double f1 = pdf[0];
        for (int i = 1; i < pdf.length; ++i) {

            double f2 = pdf[i];
            I += 0.5 * (f1 + f2) * dSigma;
            if (I >= 1.e-6) {
                return i;
            }
            f1 = f2;
        }

        return 0;
    }

    private double coverageIntervalLength(double pdf[], int i0, double dSigma) {

        double I = 0.;

        int i1 = -1;
        double f1 = pdf[i0], f2;
        for (int i = i0 + 1; i < pdf.length; ++i) {

            f2 = pdf[i];
            I += 0.5 * (f1 + f2) * dSigma;
            if (I >= p0) {
                i1 = i;
                break;
            }
            f1 = f2;
        }

        if (I < p0) { return -1.; }

        if (Math.abs(I - p0) > 1.e-3) { // should not exceed 0.1%
            System.err.println("too rough estimate, I ~ " + I);
        }

        return (i1 - i0) * dSigma;
    }

    // return a pair [start_point, length]
    private double[] coverageInterval(double pdf[], double maxSigma, double dSigma) {

        double ci[] = new double[]{-1., -1.};

        int i0 = coverageIntervalStartIndex(pdf, dSigma);

        double L0 = 10. * maxSigma;

        int nFound = 0;
        for (int i = i0; i < pdf.length; ++i) {
            double L = coverageIntervalLength(pdf, i, dSigma);
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
            ci[0] = (i0 + 1) * dSigma; //+1??
            ci[1] = L0;
        }

        return ci;
    }

    double fitM(int i, double dSigma) {

        double M = 2.;

        int it = 0, maxIt = 100;
        while (true) {

            double pdf[] = pdf(i, dSigma, (int)(M / dSigma));
            double I = dSigma * Arrays.stream(pdf).sum();

            if (Double.isNaN(I)) {
                throw new RuntimeException("I = NaN");
            }

            if (I > 0.99999) {
                break;
            }

            M = Math.round(1.5 * M);

            ++it;
            if (it > maxIt) { throw new RuntimeException(
                "cannot estimate M0 after " + maxIt + " iterations"); }
        }

        return 2. * M;
    }

    private void calculate() {

        calculatePdfParams();

        double M = fitM(0, 1.e-5);
        System.out.println("M = " + M);

        int nPts = 500_000;

        for (int i = 0; i < nSamples; ++i) {

            double dSigma = M / nPts;

            double pdf[] = pdf(i, dSigma, (int)(M / dSigma));

            // check
            double I = dSigma * Arrays.stream(pdf).sum();
            if (Math.abs(1. - I) > 1.e-4) {
                System.err.println("check failed, I ~ " + I);
                return;
            }
            double ci[] = coverageInterval(pdf, M, dSigma);
            double l = ci[0];
            double r = l + ci[1];

            if (l < 0 || r < 0) {
                System.err.println("cannot calculate CI");
                return;
            }

            System.out.printf("%2d:  CI ~ [%.4f, %.4f],  E(sigma) ~ %.4f\n",
                    i + 1, l, r, E[i]);
        }
    }


    public static void main(String args[]) { new BayesSigmaCI().calculate(); }
}
