#include <cmath>
#include <stdexcept>
#include <iostream>

// just a standalone calc. routine which does not require any libraries
class Chi2Inverse {

private:

    static constexpr double FB = 0.2316419;
    static const int NF = 5, NG = 7;
    static const double FA[], GA[], GB[], GC[];

    static double F(double z) {

        bool negZ = (z < 0.);
        if (negZ) { z *= -1.; }

        double s = 0.;
        double lambda = 1. / (1. + FB * z);
        double p = lambda;
        for (int i = 0; i < NF; ++i) {
            s += FA[i] * p;
            p *= lambda;
        }
        double c = 2.506628274631 * std::exp(0.5 * z * z);

        if (negZ) { return s / c; }
        return 1. - s / c;
    }

    static double FInv(double p) {

        if (p < 0. || p > 1) {
            throw std::invalid_argument("invalid p");
        }

        double start = -100., end = 100., mid = 0.;
        double d = 1.;
        while (d > 1.e-7) {
            mid = 0.5 * (start + end);
            double F0 = F(mid);
            if (F0 > p) { end = mid; }
            else { start = mid; }
            d = std::abs(p - F0);
        }
        return mid;
    }

public:
    // Goldstein approximation
    static double chi2Inv(double p, double n) {

        if (n < 1) { throw std::invalid_argument("invalid n"); }
        if (p < 0. || p > 1) { throw std::invalid_argument("invalid p"); }

        double s = 0.;
        double x = FInv(p);
        double tmp = 1. / std::sqrt(n);
        double v1 = 1., v2 = 1.;

        for (int i = 0; i < NG; ++i) {
            s += v1 * v2 * (GA[i] + GB[i] / n + GC[i] / (n * n));
            v1 *= x;
            v2 *= tmp;
        }
        return n * s * s * s;
    }

};

const double Chi2Inverse::FA[] = { 0.31938153, -0.35656378, 1.7814779, -1.821256, 1.3302744};
const double Chi2Inverse::GA[] = { 1.000088600,  0.471394100,  0.0001348028, -0.008553069,  0.003125580, -0.0008426812,  0.00009780499};
const double Chi2Inverse::GB[] = {-0.223736800,  0.026070830,  0.0112818600, -0.011537610,  0.005169654,  0.0025300100, -0.00145011700};
const double Chi2Inverse::GC[] = {-0.015139040, -0.008986007,  0.0227767900, -0.013232930, -0.006950356,  0.0010604380,  0.00156532600};

int main() {

    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout.precision(3);

    int dof = 10;
    std::cout << Chi2Inverse::chi2Inv(0.025, dof) << "  " << Chi2Inverse::chi2Inv(0.975, dof) << std::endl;
}
