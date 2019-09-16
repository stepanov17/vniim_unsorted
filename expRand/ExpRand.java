
import java.util.Random;

/*
 * generate random numbers from a distribution with pdf
 * f(x) ~ Math.exp( -Math.pow( Math.abs(x / (lambda * sigma)), alpha), [1]
 * lambda = Math.sqrt(gamma(1. / alpha) / gamma(3. / alpha)),
 * where sigma is a standard deviation and alpha > 0 is a parameter;
 * alpha = 2 corresponds to a normal distribution
 * alpha = +infinity corresponds to uniform distribution
 *
 * [1] Novitskiy P.V., Zograf I.A.: Ocenka pogreshnostey rezultatov izmereniy, p. 248. Energoatomizdat, Leningrad (1985) (in Russian)
 */

public class ExpRand {

    private final static double EPS = 1.e-10;

    private final static double DX = 5.e-4;

    private final static int NPERC = 99;

    private final double alpha, sigma, lambda, kf;
    private final double dx;
    private double r;
    private final Random urand;
    private final double perc[], q[];

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

    public ExpRand(double alpha, double sigma) { this(alpha, sigma, DX, null); }
    public ExpRand(double alpha, double sigma, double dx) { this(alpha, sigma, dx, null); }

    public ExpRand(double alpha, double sigma, double dx, Long seed) {

        if (alpha < EPS) { throw new RuntimeException("alpha must be positive"); }
        if (sigma < EPS) { throw new RuntimeException("sigma must be positive"); }
        if (dx < EPS) { throw new RuntimeException("dx must be positive"); }

        urand = (seed == null) ? new Random() : new Random(seed);

        this.alpha = alpha;
        this.sigma = sigma;
        this.dx = dx;

        lambda = Math.sqrt(gamma(1. / alpha) / gamma(3. / alpha));
        kf = alpha / (2. * lambda * gamma(1. / alpha) * sigma); // norming factor

        // integration range
        r = 0.;
        double f = f(r);
        while (f >= EPS) {
            r += 1.;
            f = f(r);
        }

        // percentage points
        perc = new double[NPERC];
        q = new double[NPERC];

        for (int i = NPERC; i > 0; --i) { perc[NPERC - i] = i / (NPERC + 1.); }

        for (int i = 0; i < NPERC / 2; ++i) {
            q[i] = next((i + 1.) / (NPERC + 1.), true);
            q[NPERC - i - 1] = -q[i]; // the distr. is symmetric
        }
        q[NPERC / 2] = 0.;
    }

    // pdf
    private double f(double x) {
        return kf * Math.exp( -Math.pow( Math.abs(x / (lambda * sigma)), alpha) );
    }

    private double next(double u, boolean init) {

        double x = -r, I = 0., I_prev = 0.;

        if (init) {
            if (u >= 0.5) { // symm.
                x = 0.;
                I = 0.5;
            }
        } else { // try to be a bit faster
            for (int i = 0; i < perc.length; ++i) {
                if (u >= perc[i]) {
                    x = q[i];
                    I = perc[i];
                    break;
                }
            }
        }

        for (; I < u; x += dx) {
            I_prev = I;
            I += dx * 0.125 * (f(x) + 3. * (f(x + dx / 3.) + f(x + 2. * dx / 3.)) + f(x + dx)); // 3/8 rule
        }

        // try to be a bit more accurate
        if (I > I_prev) {
            double k = (I - I_prev) / dx;
            double b = I - k * x;
            x = (u - b) / k;
        }
        return x;
    }

    // get next random value
    public double next() {
        return next(urand.nextDouble(), false);
    }


    public static void main(String args[]) {

        ExpRand rand = new ExpRand(2., 1.);

        for (int i = 0; i < 1_000_000; ++i) {
            System.out.println(i + "  " + rand.next());
        }
    }
}
