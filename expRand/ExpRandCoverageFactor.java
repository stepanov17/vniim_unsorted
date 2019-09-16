
/*
 * calculate coverage factor for a distribution with pdf
 * f(x) ~ Math.exp( -Math.pow( Math.abs(x / (lambda * sigma)), alpha), [1]
 * lambda = Math.sqrt(gamma(1. / alpha) / gamma(3. / alpha)),
 * where sigma is a standard deviation and alpha > 0 is a parameter;
 * alpha = 2 corresponds to a normal distribution
 * alpha = +infinity corresponds to uniform distribution
 *
 * [1] Novitskiy P.V., Zograf I.A.: Ocenka pogreshnostey rezultatov izmereniy, p. 248. Energoatomizdat, Leningrad (1985) (in Russian)
 */

public class ExpRandCoverageFactor {

    private final static double EPS = 1.e-10;
    private final static double DX = 1.e-4;

    private final double alpha, lambda, kf, dx;
    private double r;

    // gamma function approximation
    private double gamma(double x) {

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

    public ExpRandCoverageFactor(double alpha) { this(alpha, DX); }

    public ExpRandCoverageFactor(double alpha, double dx) {

        if (alpha < EPS) { throw new RuntimeException("alpha must be positive"); }
        if (dx < EPS) { throw new RuntimeException("dx must be positive"); }

        this.alpha = alpha;
        this.dx = dx;

        lambda = Math.sqrt(gamma(1. / alpha) / gamma(3. / alpha));
        kf = alpha / (2. * lambda * gamma(1. / alpha)); // norming factor

        // integration range
        r = 0.;
        double f = f(r);
        while (f >= EPS) {
            r += 1.;
            f = f(r);
        }
    }

    // pdf, unity variance
    private double f(double x) {
        return kf * Math.exp( -Math.pow( Math.abs(x / (lambda)), alpha) );
    }

    // get coverage factor for a given confidence level
    private double getK(double confLevel) {

        if ((confLevel < 0.7) || (confLevel > 0.999)) {
            throw new RuntimeException("please use the following confLevel range: [0.7, 0.999]"); }

        double p0 = 0.5 * confLevel; // use symmetry

        double I = 0.5, I_prev = Double.NaN;
        double x = r;

        while (I >= p0) {
            I_prev = I;
            I -= dx * 0.125 * (f(x - dx) + 3. * (f(x - 2. * dx / 3.) + f(x - dx / 3.)) + f(x)); // 3/8 rule
            x -= dx;
        }

        // try to be a bit more precise
        if ((I < p0) && (I_prev - I > EPS)) {
            double k = (I_prev - I) / dx;
            double b = I - k * x;
            double x0 = (p0 - b) / k;
            return x0;
        }

        return x;
    }

    public static void main(String args[]) {

        double P0[] = {0.75, 0.9, 0.95, 0.975, 0.99, 0.995};

        for (double p0: P0) {
            for (double alpha = 0.1; alpha < 20.0001; alpha += 0.02) {
                ExpRandCoverageFactor cf = new ExpRandCoverageFactor(alpha);
                System.out.printf("%.3f  %.3f  %.3f\n", p0, alpha, cf.getK(p0));
            }
            System.out.println("");
        }
    }
}
