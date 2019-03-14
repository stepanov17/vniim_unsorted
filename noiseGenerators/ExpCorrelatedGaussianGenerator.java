// exponentially correlated Gaussian noise
// [1] Markus Deserno, "How to generate exponentially correlated Gaussian random numbers", 2006
// (https://www.cmu.edu/biolphys/deserno/pdf/corr_gaussian_random.pdf)
//
// c = exp(-alpha), alpha > 0 - correlation between the adjacent samples
// (alpha = +infinity -> white Gaussian noise).
// noise variance should be equal to 1.

import java.util.Random;

public class ExpCorrelatedGaussianGenerator {

    private final static Random RND = new Random();
    private final double a_0, a_1;
    private double xprev;

    public ExpCorrelatedGaussianGenerator(double alpha) {
        if (alpha <= 0.) {
            throw new RuntimeException("alpha must be positive");
        }
        a_1 = Math.exp(-alpha);
        a_0 = Math.sqrt(1. - a_1 * a_1);
        xprev = 0.;
    }

    public double nextValue() {
        double x = a_0 * RND.nextGaussian() + a_1 * xprev;
        xprev = x;
        return x;
    }

    public static void main(String args[]) {
        ExpCorrelatedGaussianGenerator generator =
                new ExpCorrelatedGaussianGenerator(2.);
        for (int i = 0; i < 100; ++i) { System.out.println(generator.nextValue()); }
    }
}
