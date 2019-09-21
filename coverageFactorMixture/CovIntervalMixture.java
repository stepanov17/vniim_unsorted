
import java.util.Locale;

// calculate a shortest coverage interval for the mixture of normal distributions:
// N(0, 1) and N(m, 1), weights: w and (1 - w).
// p0 is coverage probability
public class CovIntervalMixture {

    private final double w, m;

    public CovIntervalMixture(double w, double m) {

        if (w < 0. || w > 1.) { throw new RuntimeException("invalid weight value"); }

        this.w = w;
        this.m = m;
    }

    private static class MaxIterException extends Exception {}

    private final static double A[] = {
        0.31938153, -0.35656378, 1.7814779, -1.821256, 1.3302744};
    private final static double B = 0.2316419;
    private final static double C = Math.sqrt(2. * Math.PI);

    private final static int MAX_ITER = 100_000;
    private final static double RANGE = 10.;

    // cdf approximation for normal distribution
    private double FN(double x) {

        double z = (x < 0) ? -x : x;

        double t = 1. / (1. + B * z);
        double s = 0.;
        double p = t;
        for (int i = 0; i < A.length; ++i) {
            s += A[i] * p;
            p *= t;
        }
        double v = s / (C * Math.exp(0.5 * z * z));
        return (x < 0) ? v : (1. - v);
    }

    // cdf of the mixture
    private double F(double x) {
        return w * FN(x) + (1. - w) * FN(x - m);
    }

    private double[] getCovIntervalEndPt(double x0, double dx, double p0) throws MaxIterException {

        if (p0 < 0 || p0 > 0.999) {
            throw new RuntimeException("invalid p0 value: " + p0);
        }

        double f0 = F(x0);

        double x = x0 + dx;
        double f = F(x);
        double p = f - f0, p_prev = 0.;
        int i = 0;

        while (p < p0) {
            p_prev = p;
            x += dx;
            f = F(x);
            p = f - f0;
            if (i++ > MAX_ITER) { throw new MaxIterException(); }
        }

        double res[] = new double[2];
        res[0] = x - dx;
        res[1] = p_prev;
        return res;
    }

    private double getCovIntervalEndPt(double x0, double p0) throws MaxIterException {

        double dx = 0.1;
        double x = x0, p = p0;
        while (dx > 1.e-5) {
            double tmp[] = getCovIntervalEndPt(x, dx, p);
            x = tmp[0];
            p -= tmp[1];
            dx *= 0.1;
        }
        return x;
    }

    private double[] getCovInterval(double x0, double x1, double dx, double p0) {

        double l0 = Double.POSITIVE_INFINITY;
        double v1 = Double.NaN, v2 = Double.NaN;

        for (double x = x0; x < x1 + 0.5 * dx; x += dx) {
            try {
                double y = getCovIntervalEndPt(x, p0);
                double l = y - x;
                if (l < l0) {
                    l0 = l;
                    v1 = x;
                }
            } catch (MaxIterException e) {
                break;  // just stop here
            }
        }

        if (v1 != Double.NaN) {
            v2 = v1 + l0;
        }

        double L[] = new double[2];
        L[0] = v1;
        L[1] = v2;
        return L;
    }

    private double[] getMeanAndStdev() {

        double mu = (1 - w) * m;
        double s = Math.sqrt(
                w * (1. + mu * mu) + (1. - w) * (1. + (m - mu) * (m - mu))   );
        double res[] = new double[2];
        res[0] = mu;
        res[1] = s;
        return res;
    }

    // return a pair of {cov. interval start pt, cov. interval end pt}
    private double[] getCovInterval(double p0) {

        double L[] = getCovInterval(-RANGE, RANGE, 1.e-3, p0);

        // === check ===
        double p = checkP(L[0], L[1]);
        if (Math.abs(p - p0) > 1.e-3) { throw new RuntimeException("p0 check failed, p = " + p); }

        double res[] = new double[2];
        res[0] = L[0];
        res[1] = L[1];

        return res;
    }

    // pdf
    private double f(double x) {

        double f1 = Math.exp(-0.5 * x * x) / C;
        double f2 = Math.exp(-0.5 * (x - m) * (x - m)) / C;
        return w * f1 + (1. - w) * f2;
    }

    // check coverage probability
    private double checkP(double x0, double x1) {

        if (x0 > x1) { throw new RuntimeException("check P: invalid arguments"); } 

        double I = 0., dx = 1.e-3;

        double f_prev = f(x0);

        for (double x = x0 + dx; x < x1 + 0.1 * dx; x += dx) {
            double f = f(x);
            I += 0.5 * dx * (f + f_prev);
            f_prev = f;
        }

        return I;
    }

    public static void main(String[] args) {

        double m = 3.;
        double p0 = 0.95;
        for (double w = 0.1; w < 0.901; w += 0.02) {
            CovIntervalMixture calc = new CovIntervalMixture(w, m);
            double x[] = calc.getCovInterval(p0);
            double stats[] = calc.getMeanAndStdev();
            System.out.printf(Locale.US, "%.2f  %.3f  %.3f  %.3f  %.3f  %.3f\n",
                    w, x[0], x[1], x[1] - x[0], stats[0], stats[1]);
        }
    }
}
