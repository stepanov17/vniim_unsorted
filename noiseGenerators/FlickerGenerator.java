// References:
//
// [1] Stephane Plaszczynski, Fast 1/f^alpha noise generation (https://arxiv.org/abs/astro-ph/0510081v1)
//
// [2] M. S. Keshner, 1/f noise // proceedings of the IEEE, 1982, vol. 70, issue 3
//     (https://ieeexplore.ieee.org/document/1456550)

import java.util.Random;

public class FlickerGenerator {

    // 1/f^2 generator
    private static class F2 {

        private final double a_0, a_1, b_1;
        private double xprev = 0., yprev = 0.;

        public F2(double f_min, double f_knee, double f_sample) {

            double
                r_0 = Math.PI * f_min  / f_sample,
                r_1 = Math.PI * f_knee / f_sample;

            a_0 =   (1 + r_1) / (1 + r_0);
            a_1 = - (1 - r_1) / (1 + r_0);
            b_1 =   (1 - r_0) / (1 + r_0);
        }

        public double filter(double x) {

            double y = a_0 * x + a_1 * xprev + b_1 * yprev;
            xprev = x;
            yprev = y;
            return y;
        }
    }

    // 1/f^alpha generator
    private static class F {

        private final F2 f2[]; // a set of 1/f^2 filters

        private final static Random RND = new Random();
        private final double sigma;

        public F(double alpha, double sigma, double f_min, double f_knee, double f_sample) {

            if ((alpha <= 0.) || (alpha >= 2.)) {
                System.err.println("alpha value must be in range (0 .. 2)");
                System.exit(1);
            }

            this.sigma = sigma;

            double
                    w_0 = Math.log10(2. * Math.PI * f_min),
                    w_1 = Math.log10(2. * Math.PI * f_knee);

            int N = (int)((w_1 - w_0) * 2. + Math.log10(f_sample));

            double poles[] = new double[N], zeros[] = new double[N];
            double dp = (w_1 - w_0) / N;

            poles[0] = w_0 + 0.5 * (1. - 0.5 * alpha) * dp;
            zeros[0] = poles[0] + 0.5 * alpha * dp;

            for (int i = 1; i < N; ++i) {
                poles[i] = poles[i - 1] + dp;
                zeros[i] = poles[i] + 0.5 * alpha * dp;
            }

            for (int i = 0; i < N; ++i) {
                poles[i] = Math.pow(10., poles[i]);
                zeros[i] = Math.pow(10., zeros[i]);
            }

            f2 = new F2[N];

            double c = 0.5 / Math.PI;
            for (int i = 0; i < N; ++i) {
                f2[i] = new F2(c * poles[i], c * zeros[i], f_sample);
            }
        }

        public double nextValue() {

            double v = RND.nextGaussian() * sigma;
            for (F2 f: f2) { v = f.filter(v); }
            return v;
        }
    }

    public static void main(String[] args) {

        F generator = new F(
                0.1,    // alpha
                1.,     // sigma
                1.e-5,  // f_min
                1.e-3,  // f_knee
                1.);    // f_sample

        for (int i = 0; i < 100; ++i) { System.out.println(generator.nextValue()); }
    }
}
