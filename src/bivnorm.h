// Functions distributed in this file taken from https://github.com/david-cortes/approxcdf
// Copyright 2022 David Cortes under BSD-3 license

#pragma once
#include <Rcpp.h>
#include <vector>
#include <functional>
#include <cmath>
#include <limits>
#include <algorithm>
#include <numeric>
#include <memory>
#include <cassert>


#ifndef likely
#   if defined(__GNUC__) || defined(__clang__)
#       define likely(x)   __builtin_expect(!!(x), 1)
#       define unlikely(x) __builtin_expect(!!(x), 0)
#   else
#       define likely(x)   (x)
#       define unlikely(x) (x)
#   endif
#endif
#define EPS_BLOCK 1e-20
#define LOW_RHO 1e-9
constexpr const static double HIGH_RHO = 1. - 1e-3;
/*-------------------------------------------------------*
 *  numeric constants (inline - one definition per TU)
 *-------------------------------------------------------*/
inline constexpr double inv_sqrt2_pi      = 1.0 / std::sqrt(2.0 * M_PI);
inline constexpr double neg_inv_sqrt2     = -1.0 / std::sqrt(2.0);
inline constexpr double inv_sqrt2         = 1.0 / std::sqrt(2.0);
inline constexpr double sqrt2             = std::sqrt(2.0);
inline constexpr double twoPI             = 2.0 * M_PI;
inline constexpr double fourPI            = 4.0 * M_PI;
inline constexpr double half_log_twoPI    = 0.5 * std::log(twoPI);
inline constexpr double inv_fourPI        = 1.0 / fourPI;
inline constexpr double sqrt_twoPI        = std::sqrt(2.0 * M_PI);
inline constexpr double log_sqrt_twoPI    = std::log(sqrt_twoPI);
inline constexpr double minus_inv_twoPI   = -1.0 / twoPI;
inline constexpr double inv3sqrt2pi       = 1.0 / (3.0 * std::sqrt(2.0 * M_PI));
inline constexpr double asin_one          = std::asin(1.0);
inline constexpr double inv_four_asin_one = 1.0 / (4.0 * asin_one);


/* Gauss-Legendre polynomials for quadrature rules (numerical integration) */

constexpr static const double GL5_w[] = {
    0.2369268850561891, 0.4786286704993665, 0.5688888888888889
};
constexpr static const double GL5_x[] = {
    0.9061798459386640, 0.5384693101056831, 0.0000000000000000
};

constexpr static const double GL6_w[] = {
    0.1713244923791704, 0.3607615730481386, 0.4679139345726910
};
constexpr static const double GL6_x[] = {
    0.9324695142031521, 0.6612093864662645, 0.2386191860831969
};

constexpr static const double GL8_w[] = {
    0.1012285362903763, 0.2223810344533745,
    0.3137066458778873, 0.3626837833783620
};
constexpr static const double GL8_x[] = {
    0.9602898564975363, 0.7966664774136267,
    0.5255324099163290, 0.1834346424956498
};

constexpr static const double GL12_w[] = {
    0.0471753363865118, 0.1069393259953184, 0.1600783285433462,
    0.2031674267230659, 0.2334925365383548, 0.2491470458134028
};
constexpr static const double GL12_x[] = {
    0.9815606342467192, 0.9041172563704749, 0.7699026741943047,
    0.5873179542866175, 0.3678314989981802, 0.1252334085114689
};

constexpr static const double GL16_w[] = {
    0.0271524594117541, 0.0622535239386479, 0.0951585116824928,
    0.1246289712555339, 0.1495959888165767, 0.1691565193950025,
    0.1826034150449236, 0.1894506104550685
};
constexpr static const double GL16_x[] = {
    0.9894009349916499, 0.9445750230732326, 0.8656312023878318,
    0.7554044083550030, 0.6178762444026438, 0.4580167776572274,
    0.2816035507792589, 0.0950125098376374
};

constexpr static const double GL20_w[] = {
    0.0176140071391521, 0.0406014298003869, 0.0626720483341091,
    0.0832767415767048, 0.1019301198172404, 0.1181945319615184,
    0.1316886384491766, 0.1420961093183820, 0.1491729864726037,
    0.1527533871307258
};
constexpr static const double GL20_x[] = {
    0.9931285991850949, 0.9639719272779138, 0.9122344282513259,
    0.8391169718222188, 0.7463319064601508, 0.6360536807265150,
    0.5108670019508271, 0.3737060887154195, 0.2277858511416451,
    0.0765265211334973
};

constexpr static const double GL24_w[] = {
    0.0123412297999872, 0.0285313886289337, 0.0442774388174198,
    0.0592985849154368, 0.0733464814110803, 0.0861901615319533,
    0.0976186521041139, 0.1074442701159656, 0.1155056680537256,
    0.1216704729278034, 0.1258374563468283, 0.1279381953467522
};
constexpr static const double GL24_x[] = {
    0.9951872199970213, 0.9747285559713095, 0.9382745520027328,
    0.8864155270044011, 0.8200019859739029, 0.7401241915785544,
    0.6480936519369755, 0.5454214713888396, 0.4337935076260451,
    0.3150426796961634, 0.1911188674736163, 0.0640568928626056
};

constexpr static const double GL32_w[] = {
    0.0070186100094701, 0.0162743947309057, 0.0253920653092621,
    0.0342738629130214, 0.0428358980222267, 0.0509980592623762,
    0.0586840934785355, 0.0658222227763618, 0.0723457941088485,
    0.0781938957870703, 0.0833119242269467, 0.0876520930044038,
    0.0911738786957639, 0.0938443990808046, 0.0956387200792749,
    0.0965400885147278
};
constexpr static const double GL32_x[] = {
    0.9972638618494816, 0.9856115115452684, 0.9647622555875064,
    0.9349060759377397, 0.8963211557660521, 0.8493676137325700,
    0.7944837959679424, 0.7321821187402897, 0.6630442669302152,
    0.5877157572407623, 0.5068999089322294, 0.4213512761306353,
    0.3318686022821277, 0.2392873622521371, 0.1444719615827965,
    0.0483076656877383
};

constexpr static const double GL48_w[] = {
    0.0031533460523058, 0.0073275539012763, 0.0114772345792345,
    0.0155793157229438, 0.0196161604573555, 0.0235707608393244,
    0.0274265097083569, 0.0311672278327981, 0.0347772225647704,
    0.0382413510658307, 0.0415450829434647, 0.0446745608566943,
    0.0476166584924905, 0.0503590355538545, 0.0528901894851937,
    0.0551995036999842, 0.0572772921004032, 0.0591148396983956,
    0.0607044391658939, 0.0620394231598927, 0.0631141922862540,
    0.0639242385846482, 0.0644661644359501, 0.0647376968126839
};

constexpr static const double GL48_x[] = {
    0.9987710072524261, 0.9935301722663508, 0.9841245837228269,
    0.9705915925462473, 0.9529877031604309, 0.9313866907065543,
    0.9058791367155696, 0.8765720202742479, 0.8435882616243935,
    0.8070662040294426, 0.7671590325157404, 0.7240341309238146,
    0.6778723796326639, 0.6288673967765136, 0.5772247260839727,
    0.5231609747222330, 0.4669029047509584, 0.4086864819907167,
    0.3487558862921608, 0.2873624873554556, 0.2247637903946891,
    0.1612223560688917, 0.0970046992094627, 0.0323801709628694
};

enum GLApprox {GL6=3, GL8=4, GL12=6, GL16=8, GL20=10, GL24=12, GL32=16}; 

constexpr static const double GL32_w_div4pi[] = {
    GL32_w[0] / (4. * M_PI), GL32_w[1] / (4. * M_PI), GL32_w[2] / (4. * M_PI),
    GL32_w[3] / (4. * M_PI), GL32_w[4] / (4. * M_PI), GL32_w[5] / (4. * M_PI),
    GL32_w[6] / (4. * M_PI), GL32_w[7] / (4. * M_PI), GL32_w[8] / (4. * M_PI),
    GL32_w[9] / (4. * M_PI), GL32_w[10] / (4. * M_PI), GL32_w[11] / (4. * M_PI),
    GL32_w[12] / (4. * M_PI), GL32_w[13] / (4. * M_PI), GL32_w[14] / (4. * M_PI),
    GL32_w[15] / (4. * M_PI)
};
constexpr static const double GL32_xp[] = {
    0.5 * GL32_x[0] + 0.5, 0.5 * GL32_x[1] + 0.5, 0.5 * GL32_x[2] + 0.5,
    0.5 * GL32_x[3] + 0.5, 0.5 * GL32_x[4] + 0.5, 0.5 * GL32_x[5] + 0.5,
    0.5 * GL32_x[6] + 0.5, 0.5 * GL32_x[7] + 0.5, 0.5 * GL32_x[8] + 0.5,
    0.5 * GL32_x[9] + 0.5, 0.5 * GL32_x[10] + 0.5, 0.5 * GL32_x[11] + 0.5,
    0.5 * GL32_x[12] + 0.5, 0.5 * GL32_x[13] + 0.5, 0.5 * GL32_x[14] + 0.5,
    0.5 * GL32_x[15] + 0.5
};
constexpr static const double GL32_xn[] = {
    -0.5 * GL32_x[0] + 0.5, -0.5 * GL32_x[1] + 0.5, -0.5 * GL32_x[2] + 0.5,
    -0.5 * GL32_x[3] + 0.5, -0.5 * GL32_x[4] + 0.5, -0.5 * GL32_x[5] + 0.5,
    -0.5 * GL32_x[6] + 0.5, -0.5 * GL32_x[7] + 0.5, -0.5 * GL32_x[8] + 0.5,
    -0.5 * GL32_x[9] + 0.5, -0.5 * GL32_x[10] + 0.5, -0.5 * GL32_x[11] + 0.5,
    -0.5 * GL32_x[12] + 0.5, -0.5 * GL32_x[13] + 0.5, -0.5 * GL32_x[14] + 0.5,
    -0.5 * GL32_x[15] + 0.5
};

constexpr static const double GL8_xp[] = {
    0.5 * GL8_x[0] + 0.5, 0.5 * GL8_x[1] + 0.5, 0.5 * GL8_x[2] + 0.5,
    0.5 * GL8_x[3] + 0.5
};
constexpr static const double GL8_xn[] = {
    -0.5 * GL8_x[0] + 0.5, -0.5 * GL8_x[1] + 0.5, -0.5 * GL8_x[2] + 0.5,
    -0.5 * GL8_x[3] + 0.5
};

constexpr static const double GL16_xp[] = {
    0.5 * GL16_x[0] + 0.5, 0.5 * GL16_x[1] + 0.5, 0.5 * GL16_x[2] + 0.5,
    0.5 * GL16_x[3] + 0.5, 0.5 * GL16_x[4] + 0.5, 0.5 * GL16_x[5] + 0.5,
    0.5 * GL16_x[6] + 0.5, 0.5 * GL16_x[7]
};
constexpr static const double GL16_xn[] = {
    -0.5 * GL16_x[0] + 0.5, -0.5 * GL16_x[1] + 0.5, -0.5 * GL16_x[2] + 0.5,
    -0.5 * GL16_x[3] + 0.5, -0.5 * GL16_x[4] + 0.5, -0.5 * GL16_x[5] + 0.5,
    -0.5 * GL16_x[6] + 0.5, -0.5 * GL16_x[7]
};

constexpr static const double GL48_xp[] = {
    0.5 * GL48_x[0] + 0.5, 0.5 * GL48_x[1] + 0.5, 0.5 * GL48_x[2] + 0.5,
    0.5 * GL48_x[3] + 0.5, 0.5 * GL48_x[4] + 0.5, 0.5 * GL48_x[5] + 0.5,
    0.5 * GL48_x[6] + 0.5, 0.5 * GL48_x[7] + 0.5, 0.5 * GL48_x[8] + 0.5,
    0.5 * GL48_x[9] + 0.5, 0.5 * GL48_x[10] + 0.5, 0.5 * GL48_x[11] + 0.5,
    0.5 * GL48_x[12] + 0.5, 0.5 * GL48_x[13] + 0.5, 0.5 * GL48_x[14] + 0.5,
    0.5 * GL48_x[15] + 0.5, 0.5 * GL48_x[16] + 0.5, 0.5 * GL48_x[17] + 0.5,
    0.5 * GL48_x[18] + 0.5, 0.5 * GL48_x[19] + 0.5, 0.5 * GL48_x[20] + 0.5,
    0.5 * GL48_x[21] + 0.5, 0.5 * GL48_x[22] + 0.5, 0.5 * GL48_x[23] + 0.5
};
constexpr static const double GL48_xn[] = {
    -0.5 * GL48_x[0] + 0.5, -0.5 * GL48_x[1] + 0.5, -0.5 * GL48_x[2] + 0.5,
    -0.5 * GL48_x[3] + 0.5, -0.5 * GL48_x[4] + 0.5, -0.5 * GL48_x[5] + 0.5,
    -0.5 * GL48_x[6] + 0.5, -0.5 * GL48_x[7] + 0.5, -0.5 * GL48_x[8] + 0.5,
    -0.5 * GL48_x[9] + 0.5, -0.5 * GL48_x[10] + 0.5, -0.5 * GL48_x[11] + 0.5,
    -0.5 * GL48_x[12] + 0.5, -0.5 * GL48_x[13] + 0.5, -0.5 * GL48_x[14] + 0.5,
    -0.5 * GL48_x[15] + 0.5, -0.5 * GL48_x[16] + 0.5, -0.5 * GL48_x[17] + 0.5,
    -0.5 * GL48_x[18] + 0.5, -0.5 * GL48_x[19] + 0.5, -0.5 * GL48_x[20] + 0.5,
    -0.5 * GL48_x[21] + 0.5, -0.5 * GL48_x[22] + 0.5, -0.5 * GL48_x[23] + 0.5
};
#ifndef _OPENMP
constexpr static const double GL5_w_div4pi[] = {
    GL5_w[0] / (4. * M_PI), GL5_w[1] / (4. * M_PI), GL5_w[2] / (4. * M_PI)
};
constexpr static const double GL5_xp[] = {
    0.5 * GL5_x[0] + 0.5, 0.5 * GL5_x[1] + 0.5, 0.5 * GL5_x[2] + 0.5
};
constexpr static const double GL5_xn[] = {
    -0.5 * GL5_x[0] + 0.5, -0.5 * GL5_x[1] + 0.5, -0.5 * GL5_x[2] + 0.5
};
#else
constexpr static const double GL8_w_div4pi[] = {
    GL8_w[0] / (4. * M_PI), GL8_w[1] / (4. * M_PI),
    GL8_w[2] / (4. * M_PI), GL8_w[3] / (4. * M_PI)
};
#endif
double norm_pdf_1d(double x)
{
    if (std::isinf(x)) return 0.;
    if (likely(x < 5. && x > -5.)) {
        return inv_sqrt2_pi * std::exp(-0.5 * x * x);
    }
    else {
        double x1 = std::floor(x * 0x1.0p16 + 0.5) * 0x1.0p-16;
        double x2 = x - x1;
        return inv_sqrt2_pi * (std::exp(-0.5 * x1 * x1) * std::exp((-0.5 * x2 - x1) * x2));
    }
}

/* Adapted from cephes' 'ndtr' */
double norm_cdf_1d(double x)
{
    /* This is technically correct but very imprecise:  
    return 0.5 * std::erfc(neg_inv_sqrt2 * x);
    */
    if (std::isinf(x)) {
        return (x >= 0.)? 1. : 0.;
    }
    double x_, y, z;
    x_ = x * inv_sqrt2;
    z = std::fabs(x_);

    /* if( z < SQRTH ) */
    // if (z < 1.) {
    if (z < inv_sqrt2) {
        y = .5 + .5*std::erf(x_);
    }
    else {
        y = .5*std::erfc(z);
        if(x_ > 0.) {
            y = 1. - y;
        }
    }
    return y;
}

double norm_lcdf_1d(double x)
{
    /* This is technically correct but very imprecise:
    return 0.5 * std::erfc(inv_sqrt2 * x);
    */
    return norm_cdf_1d(-x);
}

/* Adapted from SciPy:
   https://github.com/scipy/scipy/blob/main/scipy/stats/_continuous_distns.py */
double norm_logpdf_1d(double x)
{
    if (std::isinf(x)) return -std::numeric_limits<double>::infinity();
    return -0.5 * (x*x) - log_sqrt_twoPI;
}

/* Adapted from SciPy:
   https://github.com/scipy/scipy/blob/8a64c938ddf1ae4c02a08d2c5e38daeb8d061d38/scipy/special/cephes/ndtr.c */
double norm_logcdf_1d(double a)
{
    if (std::isinf(a)) {
        return (a >= 0.)? 0. : -std::numeric_limits<double>::infinity();
    }
    const double a_sq = a * a;
    double log_LHS;              /* we compute the left hand side of the approx (LHS) in one shot */
    double last_total = 0;       /* variable used to check for convergence */
    double right_hand_side = 1;  /* includes first term from the RHS summation */
    double numerator = 1;        /* numerator for RHS summand */
    double denom_factor = 1;     /* use reciprocal for denominator to avoid division */
    double denom_cons = 1./a_sq; /* the precomputed division we use to adjust the denominator */
    long sign = 1;
    long i = 0;

    if (a > 6.) {
        return -norm_cdf_1d(-a);        /* log(1+x) \approx x */
    }
    if (a > -20.) {
        return std::log(norm_cdf_1d(a));
    }
    log_LHS = -0.5*a_sq - std::log(-a) - half_log_twoPI;

    while (std::fabs(last_total - right_hand_side) > std::numeric_limits<double>::epsilon()) {
        i++;
        last_total = right_hand_side;
        sign = -sign;
        denom_factor *= denom_cons;
        numerator *= 2 * i - 1;
        right_hand_side = std::fma(sign*numerator, denom_factor, right_hand_side);
    }
    
    return log_LHS + std::log(right_hand_side);
}

/* Bivariate normal CDF.
   Algorithm:  Drezner (1978) 5-point GL with the p-split
               refined by Drezner & Wesolowsky (1990);
               parameter cut-off as in West (2004).
   Expected accuracy ~ 1e-6.
*/
double norm_lcdf_2d_fast(double x1, double x2, double rho)
{
    double x12 = 0.5 * (x1*x1 + x2*x2);
    
    double out = 0;
    double r1, x3;
    if (std::fabs(rho) >= 0.7) {
        double r2 = 1. - rho*rho;
        double r3 = std::sqrt(r2);
        if (rho < 0) {
            x2 = -x2;
        }
        x3 = x1*x2;
        double x7 = std::exp(-0.5 * x3);
        if (r2) {
            double x6 = std::fabs(x1 - x2);
            double x5 = 0.5 * x6*x6;
            x6 /= r3;
            double aa = 0.5 - 0.125*x3;
            double ab = 3. - 2.*aa*x5;
            out = inv3sqrt2pi * (
                x6 * ab * norm_lcdf_1d(x6) -
                std::exp(-x5/r2) * std::fma(aa, r2, ab) * inv_sqrt2_pi
            );
            double rr;
            double nr1, nrr, nr2;
            #ifndef _OPENMP
            #pragma GCC unroll 2
            for (int ix = 0; ix < 2; ix++) {
                r1 = r3 * GL5_xp[ix];
                rr = r1*r1;
                r2 = std::sqrt(1. - rr);

                nr1 = r3 * GL5_xn[ix];
                nrr = nr1*nr1;
                nr2 = std::sqrt(1. - nrr);
                
                out -= GL5_w_div4pi[ix] * (
                    std::exp(-x5/rr) * (std::exp(-x3/(1. + r2))/r2/x7 - 1.- aa*rr) +
                    std::exp(-x5/nrr) * (std::exp(-x3/(1. + nr2))/nr2/x7 - 1.- aa*nrr)
                );
            }
            r1 = r3 * GL5_xp[2];
            rr = r1*r1;
            r2 = std::sqrt(1. - rr);
            out -= GL5_w_div4pi[2] * std::exp(-x5/rr) * (std::exp(-x3/(1. + r2))/r2/x7 - 1.- aa*rr);
            #else
            #ifndef _MSC_VER
            #pragma omp simd
            #endif
            for (int ix = 0; ix < 4; ix++) {
                r1 = r3 * GL8_xp[ix];
                rr = r1*r1;
                r2 = std::sqrt(1. - rr);

                nr1 = r3 * GL8_xn[ix];
                nrr = nr1*nr1;
                nr2 = std::sqrt(1. - nrr);
                
                out -= GL8_w_div4pi[ix] * (
                    std::exp(-x5/rr) * (std::exp(-x3/(1. + r2))/r2/x7 - 1.- aa*rr) +
                    std::exp(-x5/nrr) * (std::exp(-x3/(1. + nr2))/nr2/x7 - 1.- aa*nrr)
                );
            }
            #endif
        }
        if (rho > 0) {
            out = std::fma(out, r3*x7, norm_lcdf_1d(std::fmax(x1, x2)));
        }
        else {
            out = std::fmax(0., norm_lcdf_1d(x1) - norm_lcdf_1d(x2)) - out*r3*x7;
        }
        return out;
    }
    else {
        x3 = x1*x2;
        double rr2;
        double nr1, nrr2;
        #ifndef _OPENMP
        #pragma GCC unroll 2
        for (int ix = 0; ix < 2; ix++) {
            r1 = rho * GL5_xp[ix];
            rr2 = 1. - r1*r1;

            nr1 = rho * GL5_xn[ix];
            nrr2 = 1. - nr1*nr1;

            out += GL5_w_div4pi[ix] * (
                std::exp((r1*x3 - x12) / rr2) / std::sqrt(rr2) +
                std::exp((nr1*x3 - x12) / nrr2) / std::sqrt(nrr2)
            );
        }
        r1 = rho * GL5_xp[2];
        rr2 = 1. - r1*r1;
        out += GL5_w_div4pi[2] * std::exp((r1*x3 - x12) / rr2) / std::sqrt(rr2);
        #else
        #ifndef _MSC_VER
        #pragma omp simd
        #endif
        for (int ix = 0; ix < 4; ix++) {
            r1 = rho * GL8_xp[ix];
            rr2 = 1. - r1*r1;

            nr1 = rho * GL8_xn[ix];
            nrr2 = 1. - nr1*nr1;

            out += GL8_w_div4pi[ix] * (
                std::exp((r1*x3 - x12) / rr2) / std::sqrt(rr2) +
                std::exp((nr1*x3 - x12) / nrr2) / std::sqrt(nrr2)
            );
        }
        #endif
        return std::fma(out, rho, norm_lcdf_1d(x1) * norm_lcdf_1d(x2));
    }
}

double norm_cdf_2d_fast(double x1, double x2, double rho)
{
    return norm_lcdf_2d_fast(-x1, -x2, rho);
}

/* Tsay, Wen-Jen, and Peng-Hsuan Ke.
   "A simple approximation for the bivariate normal integral."
   Communications in Statistics-Simulation and Computation (2021): 1-14. */
constexpr const static double c1 = -1.0950081470333;
constexpr const static double c2 = -0.75651138383854;
double norm_cdf_2d_vfast(double x1, double x2, double rho)
{
    if (std::fabs(rho) <= std::numeric_limits<double>::epsilon()) {
        return norm_cdf_1d(x1) * norm_cdf_1d(x2);
    }

    double denom = std::sqrt(1 - rho * rho);
    double a = -rho / denom;
    double b = x1 / denom;
    double aq_plus_b = a*x2 + b;

    if (a > 0) {
        if (aq_plus_b >= 0) {
            double aa = a * a;
            double a_sq_c1 = aa*c1;
            double a_sq_c2 = aa*c2;
            double sqrt2b = sqrt2*b;
            double sqrt2x2 = sqrt2*x2;
            double sqrt_recpr_a_sq_c2 = std::sqrt(1. - a_sq_c2);
            double twicea_sqrt_recpr_a_sq_c2 = 2.*a*sqrt_recpr_a_sq_c2;
            double temp = 1. / (4. * sqrt_recpr_a_sq_c2);
            double t1 = a_sq_c1*c1 + 2.*b*b*c2;
            double t2 = 2.*sqrt2b*c1;
            double t3 = 4. - 4.*a_sq_c2;

            return
                0.5 * (std::erf(x2 / sqrt2) + std::erf(b / (sqrt2*a))) +
                temp
                    * std::exp((t1 - t2) / t3)
                    * (1. - std::erf((sqrt2b - a_sq_c1) / twicea_sqrt_recpr_a_sq_c2)) -
                temp
                    * std::exp((t1 + t2) / t3)
                    * (
                        std::erf((sqrt2x2 - sqrt2x2*a_sq_c2 - sqrt2b*a*c2 - a*c1) / (2.*sqrt_recpr_a_sq_c2)) +
                        std::erf((a_sq_c1 + sqrt2b) / twicea_sqrt_recpr_a_sq_c2)
                    );

        }
        else {
            double sqrt2b = sqrt2*b;
            double sqrt2x2 = sqrt2*x2;
            double a_sq_c2 = a*a*c2;
            double recpr_a_sq_c2 = 1. - a_sq_c2;
            double sqrt_recpr_a_sq_c2 = std::sqrt(recpr_a_sq_c2);
            double a_c1 = a*c1;

            return
                (1. / (4. * sqrt_recpr_a_sq_c2)) *
                std::exp((a_c1*a_c1 - 2.*sqrt2b*c1 + 2*b*b*c2) / (4.*recpr_a_sq_c2)) *
                (1. + std::erf((sqrt2x2 - sqrt2x2*a_sq_c2 - sqrt2b*a*c2 + a_c1) / (2.*sqrt_recpr_a_sq_c2)));
        }
    }
    else {
        if (aq_plus_b >= 0) {
            double sqrt2b = sqrt2*b;
            double a_sq_c2 = a*a*c2;
            double recpr_a_sq_c2 = 1. - a_sq_c2;
            double sqrt_recpr_a_sq_c2 = std::sqrt(recpr_a_sq_c2);
            double a_c1 = a*c1;
            double sqrt2_x2 = sqrt2*x2;

            return
                0.5 + 0.5 * std::erf(x2 / sqrt2) -
                (1. / (4. * sqrt_recpr_a_sq_c2)) *
                std::exp((a_c1*a_c1 + 2.*sqrt2b*c1 + 2.*b*b*c2) / (4.*recpr_a_sq_c2)) *
                (1. + std::erf((sqrt2_x2 - sqrt2_x2*a_sq_c2 - sqrt2b*a*c2 - a_c1) / (2.*sqrt_recpr_a_sq_c2)));
        }
        else {
            double sqrt2a = sqrt2*a;
            double sqrt2b = sqrt2*b;
            double a_sq_c2 = a*a*c2;
            double recpr_a_sq_c2 = 1. - a_sq_c2;
            double sqrt_recpr_a_sq_c2 = std::sqrt(recpr_a_sq_c2);
            double a_c1 = a*c1;
            double temp = 1. / (4. * sqrt_recpr_a_sq_c2);
            double t1 = a_c1*a_c1 + 2.*b*b*c2;
            double t2 = 2.*sqrt2b*c1;
            double t3 = 4.*recpr_a_sq_c2;
            double sqrt2_x2 = sqrt2*x2;

            return
                0.5 - 0.5 * std::erf(b / sqrt2a) -
                temp
                    * std::exp((t1 + t2) / t3)
                    * (1. - std::erf((sqrt2b + a*a_c1) / (2.*a*sqrt_recpr_a_sq_c2))) +
                temp
                    * std::exp((t1 - t2) / t3)
                    * (
                        std::erf((sqrt2_x2 - sqrt2_x2*a_sq_c2 - sqrt2b*a*c2 + a_c1) / (2.*sqrt_recpr_a_sq_c2)) +
                        std::erf((sqrt2b - a*a_c1) / (2.*a*sqrt_recpr_a_sq_c2))
                    );
        }
    }
} 

[[gnu::flatten]]
double norm_lcdf_2d(double x1, double x2, double rho)
{
    double abs_rho = std::fabs(rho);

    double out = 0;
    double hk;
    if (abs_rho < 0.925) {
        if (abs_rho > std::numeric_limits<double>::epsilon()) {
            hk = x1 * x2;
            double hs = 0.5 * (x1*x1 + x2*x2);
            double asr = std::asin(rho);
            double asr_half = 0.5 * asr;
            double sn1, sn2;

            int gl_dim;
            #ifndef _OPENMP
            const double* GL_wtable;
            const double* GL_xtable;
            #endif
            if (abs_rho < 0.3) {
                gl_dim = (int)GL6;
                #ifndef _OPENMP
                GL_wtable = GL6_w;
                GL_xtable = GL6_x;
                #endif
            }
            else if (abs_rho < 0.5) {
                gl_dim = (int)GL12;
                #ifndef _OPENMP
                GL_wtable = GL12_w;
                GL_xtable = GL12_x;
                #endif
            }
            else {
                gl_dim = (int)GL20;
                #ifndef _OPENMP
                GL_wtable = GL20_w;
                GL_xtable = GL20_x;
                #endif
            }
            
            #ifndef _OPENMP
            #pragma GCC unroll 3
            for (int ix = 0; ix < gl_dim; ix++) {
                sn1 = std::sin(asr_half * (1. + GL_xtable[ix]));
                sn2 = std::sin(asr_half * (1. - GL_xtable[ix]));
                out += GL_wtable[ix] * (
                    std::exp(std::fma(sn1, hk, -hs) / std::fma(-sn1, sn1, 1.)) +
                    std::exp(std::fma(sn2, hk, -hs) / std::fma(-sn2, sn2, 1.))
                );
            }
            #else
            switch (gl_dim) {
                case GL6: {
                    #ifndef _MSC_VER
                    #pragma omp simd
                    #endif
                    for (int ix = 0; ix < 4; ix++) {
                        sn1 = std::sin(asr_half * (1. + GL8_x[ix]));
                        sn2 = std::sin(asr_half * (1. - GL8_x[ix]));
                        out += GL8_w[ix] * (
                            std::exp(std::fma(sn1, hk, -hs) / std::fma(-sn1, sn1, 1.)) +
                            std::exp(std::fma(sn2, hk, -hs) / std::fma(-sn2, sn2, 1.))
                        );
                    }
                    break;
                }
                case GL12: {
                    #ifndef _MSC_VER
                    #pragma omp simd
                    #endif
                    for (int ix = 0; ix < 8; ix++) {
                        sn1 = std::sin(asr_half * (1. + GL16_x[ix]));
                        sn2 = std::sin(asr_half * (1. - GL16_x[ix]));
                        out += GL16_w[ix] * (
                            std::exp(std::fma(sn1, hk, -hs) / std::fma(-sn1, sn1, 1.)) +
                            std::exp(std::fma(sn2, hk, -hs) / std::fma(-sn2, sn2, 1.))
                        );
                    }
                    break;
                }
                case GL20: {
                    #ifndef _MSC_VER
                    #pragma omp simd
                    #endif
                    for (int ix = 0; ix < 12; ix++) {
                        sn1 = std::sin(asr_half * (1. + GL24_x[ix]));
                        sn2 = std::sin(asr_half * (1. - GL24_x[ix]));
                        out += GL24_w[ix] * (
                            std::exp(std::fma(sn1, hk, -hs) / std::fma(-sn1, sn1, 1.)) +
                            std::exp(std::fma(sn2, hk, -hs) / std::fma(-sn2, sn2, 1.))
                        );
                    }
                    break;
                }
            }
            #endif
            out *= asr / fourPI;
        }
        out = std::fma(norm_lcdf_1d(x1), norm_lcdf_1d(x2), out);
    }
    else {

        if (rho < 0) {
            x2 = -x2;
        }
        if (abs_rho < 1) {
            hk = x1 * x2;
            double as = std::fma(-rho, rho, 1.);
            double a = std::sqrt(as);
            double b;
            double bs = (x1 - x2) * (x1 - x2);
            double c = std::fma(-0.125, hk, 0.5);
            double d = std::fma(-0.0625, hk, 0.75);
            double asr = -0.5 * (hk + bs / as);
            double rfdbs = std::fma(-d, bs, 5.)*(1./15.);
            if (asr > -100.) {
                out = a * std::exp(asr) * (1. - c * (bs - as) * rfdbs + 0.2*c*d*as*as);
            }
            if (hk > -100.) {
                b = std::sqrt(bs);
                out -= std::exp(-0.5 * hk) * sqrt_twoPI * norm_lcdf_1d(b / a) * b * (1. - c*bs*rfdbs);
            }
            a *= 0.5;
            double xs;
            double rs;
            double temp;

            #ifndef _OPENMP
            #pragma GCC unroll 10
            for (int ix = 0; ix < 10; ix++) {
                temp = a * (1. + GL20_x[ix]);
                xs = temp * temp;
                rs = std::sqrt(1. - xs);
                asr = -0.5 * (hk + bs / xs);
                if (asr > -100.) {
                    out += a * GL20_w[ix] * std::exp(asr) *
                           (std::exp(-hk*xs/(2.*(1.+rs)*(1.+rs)))/rs - (1. + c*xs*std::fma(d, xs, 1.)));
                }

                temp = a * (1. - GL20_x[ix]);
                xs = temp * temp;
                rs = std::sqrt(1. - xs);
                asr = -0.5 * (hk + bs / xs);
                if (asr > -100.) {
                    out += a * GL20_w[ix] * std::exp(asr) *
                           (std::exp(-hk*xs/(2.*(1.+rs)*(1.+rs)))/rs - (1. + c*xs*std::fma(d, xs, 1.)));
                }
            }
            #else
            #ifndef _MSC_VER
            #pragma omp simd
            #endif
            for (int ix = 0; ix < 12; ix++) {
                temp = a * (1. + GL24_x[ix]);
                xs = temp * temp;
                rs = std::sqrt(1. - xs);
                asr = -0.5 * (hk + bs / xs);
                out += a * GL24_w[ix] * std::exp(asr) *
                       (std::exp(-hk*xs/(2.*(1.+rs)*(1.+rs)))/rs - (1. + c*xs*std::fma(d, xs, 1.)));

                temp = a * (1. - GL24_x[ix]);
                xs = temp * temp;
                rs = std::sqrt(1. - xs);
                asr = -0.5 * (hk + bs / xs);
                out += a * GL24_w[ix] * std::exp(asr) *
                       (std::exp(-hk*xs/(2.*(1.+rs)*(1.+rs)))/rs - (1. + c*xs*std::fma(d, xs, 1.)));
            }
            #endif
            out *= minus_inv_twoPI;
        }
        if (rho > 0) {
            out += norm_lcdf_1d(std::fmax(x1, x2));
        }
        else {
            out = -out;
            if (x2 > x1) {
                if (x1 < 0) {
                    out += norm_cdf_1d(x2) - norm_cdf_1d(x1);
                }
                else {
                    out += norm_lcdf_1d(x1) - norm_lcdf_1d(x2);
                }
            }
        }
    }
    return out;
}

double norm_cdf_2d(double x1, double x2, double rho)
{
    return norm_lcdf_2d(-x1, -x2, rho);
}

// [[Rcpp::export]]
double pbvn_tsay(double h, double k, double rho) {
    return norm_cdf_2d_vfast(h, k, rho);
}

// [[Rcpp::export]]
double pbvn_tvpack(double h, double k, double rho) {
    return norm_cdf_2d(h, k, rho);
}

// [[Rcpp::export]]
double pbvn_drezner(double h, double k, double rho) {
    return norm_cdf_2d_fast(h, k, rho);
}
