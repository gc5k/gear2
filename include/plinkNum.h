#ifndef PLINKNUM_H
#define PLINKNUM_H

#include "time.h"
#include <bits/stdc++.h>
#include <math.h>
using namespace std;

namespace plinknum {

//***********PLINK_DEFINE*********
#define S_CAST(type, val) (static_cast<type>(val))
#define MINV(aa, bb) (((bb) < (aa))? (bb) : (aa))
//***********PLINK_DEFINE*********

//***********PLINK_CONST*********
static const double kLogMinValue = -708.0;
static const double kLanczosSumNumer[6] = {8706.3495925490091, 8523.650341121874, 3338.0292194764235, 653.64249942940087, 63.999518449381870, 2.5066282746310063};
static const double kLanczosSumDenom[6] = {0, 24, 50, 35, 10, 1};
static const double kLanczosSumExpgNumer[6] = {32.812445410297834, 32.123889414443320, 12.580347294552161, 2.4634444783532414, 0.2412010548258800, 0.0094469677045392};
static const double kLanczosG = 5.581;
static const double kEpsilon = 0.000000000931322574615478515625;
static const double kLnSqrtPi = 0.5723649429247001;
// 2^{-21}, must be >= sqrt(kSmallEpsilon)
static const double kBigEpsilon = 0.000000476837158203125;
// they can't live in plink2_cmdline any more.
static const double kE = 2.7182818284590452;
static const double kPi = 3.1415926535897932;
static const double kSqrt2 = 1.4142135623730951;
static const double kRecipE = 0.36787944117144233;
static const double kRecip2m53 = 0.00000000000000011102230246251565404236316680908203125;
static const double kSqrtPi = 1.7724538509055159;

static const double kLn2 = 0.6931471805599453;
static const double kLentzFpmin = 1.0e-30;
// ~6 digits of precision is appropriate for p-value computations
static const double kContinuedFractionEpsilon = 3.0e-7;
static const double kExactTestBias = 0.00000000000000000000000010339757656912845935892608650874535669572651386260986328125;
static const double kSmallEpsilon = 0.00000000000005684341886080801486968994140625;

static const double kTemmeC0[7] = {-0.333333333, 0.0833333333, -0.0148148148, 0.00115740741, 0.000352733686, -0.000178755144, 0.391926318e-4};
static const double kTemmeC1[5] = {-0.00185185185, -0.00347222222, 0.00264550265, -0.000990226337, 0.000205761317};
static const double kTemmeC2[3] = {0.00413359788, -0.00268132716, 0.000771604938};

static const double kSmallHalfRecips[30] = {
  1.0 / 0.5,
  1.0 / 1.5,
  1.0 / 2.5,
  1.0 / 3.5,
  1.0 / 4.5,
  1.0 / 5.5,
  1.0 / 6.5,
  1.0 / 7.5,
  1.0 / 8.5,
  1.0 / 9.5,
  1.0 / 10.5,
  1.0 / 11.5,
  1.0 / 12.5,
  1.0 / 13.5,
  1.0 / 14.5,
  1.0 / 15.5,
  1.0 / 16.5,
  1.0 / 17.5,
  1.0 / 18.5,
  1.0 / 19.5,
  1.0 / 20.5,
  1.0 / 21.5,
  1.0 / 22.5,
  1.0 / 23.5,
  1.0 / 24.5,
  1.0 / 25.5,
  1.0 / 26.5,
  1.0 / 27.5,
  1.0 / 28.5,
  1.0 / 29.5
};

static const double kFactorialRecips[34] = {
  1.0,
  1.0,
  0.5,
  0.16666666666666666,
  0.041666666666666664,
  0.008333333333333333,
  0.001388888888888889,
  0.0001984126984126984,
  2.48015873015873e-05,
  2.7557319223985893e-06,
  2.755731922398589e-07,
  2.505210838544172e-08,
  2.08767569878681e-09,
  1.6059043836821613e-10,
  1.1470745597729725e-11,
  7.647163731819816e-13,
  4.779477332387385e-14,
  2.8114572543455206e-15,
  1.5619206968586225e-16,
  8.22063524662433e-18,
  4.110317623312165e-19,
  1.9572941063391263e-20,
  8.896791392450574e-22,
  3.8681701706306835e-23,
  1.6117375710961184e-24,
  6.446950284384474e-26,
  2.4795962632247972e-27,
  9.183689863795546e-29,
  3.279889237069838e-30,
  1.1309962886447718e-31,
  3.7699876288159054e-33,
  1.2161250415535181e-34,
  3.800390754854744e-36,
  1.151633562077195e-37
};

static const double kSmallRecips[30] = {
  0.0,  // could make this nan, though that's annoying for C++03
  1.0 / 1,
  1.0 / 2,
  1.0 / 3,
  1.0 / 4,
  1.0 / 5,
  1.0 / 6,
  1.0 / 7,
  1.0 / 8,
  1.0 / 9,
  1.0 / 10,
  1.0 / 11,
  1.0 / 12,
  1.0 / 13,
  1.0 / 14,
  1.0 / 15,
  1.0 / 16,
  1.0 / 17,
  1.0 / 18,
  1.0 / 19,
  1.0 / 20,
  1.0 / 21,
  1.0 / 22,
  1.0 / 23,
  1.0 / 24,
  1.0 / 25,
  1.0 / 26,
  1.0 / 27,
  1.0 / 28,
  1.0 / 29,
  
};
//***********PLINK_CONST*********



	//***********PLINK_FUNCTION*********
	double u31tod(uint32_t uii) {
		const int32_t ii = uii;
		assert(ii >= 0);
		return S_CAST(double, ii);
	}

  double lanczos_sum_expg_scaled_imp(double zz, double *s2_ptr)
  {
    double s1;
    double s2;
    if (zz <= 1)
    {
      s1 = kLanczosSumExpgNumer[5];
      s2 = kLanczosSumDenom[5];
      for (int32_t ii = 4; ii >= 0; --ii)
      {
        s1 *= zz;
        s2 *= zz;
        s1 += kLanczosSumExpgNumer[S_CAST(uint32_t, ii)];
        s2 += kLanczosSumDenom[S_CAST(uint32_t, ii)];
      }
    }
    else
    {
      zz = 1 / zz;
      s1 = kLanczosSumExpgNumer[0];
      s2 = kLanczosSumDenom[0];
      for (uint32_t uii = 1; uii != 6; ++uii)
      {
        s1 *= zz;
        s2 *= zz;
        s1 += kLanczosSumExpgNumer[uii];
        s2 += kLanczosSumDenom[uii];
      }
    }
    *s2_ptr = s2;
    return s1;
  }

double ibeta_series_ln(double aa, double bb, double xx, uint32_t inv) {
  // BPSER in DiDonato and Morris.
  // ibeta(a, b, x) ~= gamma(a+b) / (gamma(a) * gamma(b))
  //                 * x^a / a
  //                 * (1 + a\sum_j=1^N [(1-b)(2-b)...(j-b)x^j / (j!(a+j))])
  // Currently always called with x <= 0.7, and either b*x <= 0.7 or b<40, so
  // guaranteed to converge at reasonable speed.

  // normalized always true

  const double cc = aa + bb;

  const double agh = aa + kLanczosG - 0.5;
  const double bgh = bb + kLanczosG - 0.5;
  const double cgh = cc + kLanczosG - 0.5;
  double numer_a;
  double denom_a = lanczos_sum_expg_scaled_imp(aa, &numer_a);
  double numer_b;
  double denom_b = lanczos_sum_expg_scaled_imp(bb, &numer_b);
  double denom_c;
  double numer_c = lanczos_sum_expg_scaled_imp(cc, &denom_c);
  double result = (numer_a * numer_b * numer_c) / (denom_a * denom_b * denom_c);
  double l1 = log(cgh / bgh) * (bb - 0.5);
  double l2 = log(xx * cgh / agh) * aa;
  double result_ln = log(result * result * agh) * 0.5 + l1 + l2 - 0.5;

  double series_sum = 0.0;
  double term = 1.0;
  double apn = aa;
  double poch = 1.0 - bb;
  double nn = 1.0;
  while (1) {
    double rr = term / apn;
    series_sum += rr;
    if (fabs(rr) <= fabs(series_sum * kEpsilon)) {
      result_ln += log(series_sum);
      if (!inv) {
        return result_ln;
      }
      // assume for now that we wouldn't be inverting into the tiny result_ln
      // case.
      return log1p(-exp(result_ln));
    }
    apn += 1.0;
    term *= poch * xx / nn;
    nn += 1.0;
    poch += 1.0;
  }
}

double neg_powm1_imp_ln(double xx, double yy) {
  const double ll = yy * log(xx);
  if (ll > -1.0) {
    // For tiny |l|, we may lose all precision when exp(l) evaluates to 1.
    return log(-expm1(ll));
  }
  // For larger |l|, it's the log step that may throw away precision.
  if (ll < kLogMinValue) {
    return 0.0;
  }
  return log1p(-exp(ll));
}

double log1pmx(double xx) {
  // log(1+x) - x
  // assumes abs(xx) < 0.95
  const double aa = fabs(xx);
  if (aa < (kBigEpsilon / kSqrt2)) { // 2^{-21.5}
    return -xx * xx * 0.5;
  }
  double kk = 1.0;  // skip first term of usual log(1+x) series
  const double m_mult = -xx;
  double m_prod = xx;
  double total = 0.0;
  double rr;
  do {
    m_prod *= m_mult;
    kk += 1.0;
    rr = m_prod / kk;
    total += rr;
    // todo: tune these epsilons, but let's wait until we know all of the
    // callers of these functions
  } while (fabs(rr) > (kBigEpsilon * kBigEpsilon));
  return total;
}

double lanczos_sum_expg_scaled_recip(double zz) {
  double s2;
  double s1 = lanczos_sum_expg_scaled_imp(zz, &s2);
  return s2 / s1;
}

// compute -log((z^a)(e^{-z})/tgamma(a))
double regularized_gamma_prefix_ln(double aa, double zz) {
  // assumes a == 0.5 if a < 1.  assumes z > 0.
  // we are fine with float-level precision, so lanczos_n=6, kLanczosG=5.581
  if (aa < 1) {
    return -zz + 0.5 * log(zz) - kLnSqrtPi;
  }
  const double agh = aa + kLanczosG - 0.5;
  const double agh_recip = 1.0 / agh;
  const double dd = ((zz - aa) - (kLanczosG - 0.5)) * agh_recip;
  double prefix_ln;
  if ((fabs(dd * dd * aa) <= 100) && (aa > 150)) {
    // abs(dd) < sqrt(2/3) < 0.95
    prefix_ln = aa * log1pmx(dd) + zz * (0.5 - kLanczosG) * agh_recip;
  } else {
    prefix_ln = (aa - zz) + aa * log(zz * agh_recip);
  }
  const double scaled_recip = lanczos_sum_expg_scaled_recip(aa);
  return prefix_ln + 0.5 * log(agh * kRecipE * scaled_recip * scaled_recip);
}

double lanczos_sum(double zz) {
  double s1;
  double s2;
  if (zz <= 1) {
    s1 = kLanczosSumNumer[5];
    s2 = kLanczosSumDenom[5];
    for (int32_t ii = 4; ii >= 0; --ii) {
      s1 *= zz;
      s2 *= zz;
      s1 += kLanczosSumNumer[S_CAST(uint32_t, ii)];
      s2 += kLanczosSumDenom[S_CAST(uint32_t, ii)];
    }
  } else {
    zz = 1 / zz;
    s1 = kLanczosSumNumer[0];
    s2 = kLanczosSumDenom[0];
    for (uint32_t uii = 1; uii != 6; ++uii) {
      s1 *= zz;
      s2 *= zz;
      s1 += kLanczosSumNumer[uii];
      s2 += kLanczosSumDenom[uii];
    }
  }
  return s1 / s2;
}

double tgamma_delta_ratio(double zz, double delta) {
  // gamma(z) / gamma(z + delta)
  // We only call this with delta = 0.5 or 1, so no overflow issues.
  // zz >= 15 for now.
  if (delta == 1.0) {
    // Trivial case.
    return 1.0 / zz;
  }
  assert(delta == 0.5);
  // Can skip z < epsilon and z + delta == z branches for now.
  double zgh = zz + kLanczosG - 0.5;
  // Also skip fabs(delta) >= 10 branch for now.
  double result = exp((0.5 - zz) * log1p(delta / zgh));
  result *= lanczos_sum(zz) / lanczos_sum(zz + delta);
  // exploit forced delta == 0.5
  result *= sqrt(kE / (zgh + delta));
  return result;
}

double erfc_fast2(double zz, double* tau_ln_plus_z2_ptr) {
  const double tt = 1.0 / (1.0 + 0.5 * zz);
  *tau_ln_plus_z2_ptr = ((((((((0.17087277 * tt - 0.82215223) * tt + 1.48851587) * tt - 1.13520398) * tt + 0.27886807) * tt - 0.18628806) * tt + 0.09678418) * tt + 0.37409196) * tt + 1.00002368) * tt - 1.26551223;
  return tt;
}

double finite_half_gamma_q2_ln(uint32_t a_minus_half, double xx) {
  // a is in {0.5, 1.5, ..., 29.5}; max(0.2, a-1) < x < 1e10
  const double sqrt_x = sqrt(xx);
  double tau_ln_plus_x;
  double tt = erfc_fast2(sqrt_x, &tau_ln_plus_x);
  if (!a_minus_half) {
    return log(tt) + tau_ln_plus_x - xx;
  }
  // pre-multiply by e^x to avoid underflow
  double term = sqrt_x * (2.0 / kSqrtPi);
  double sum = term;
  for (uint32_t n_minus_half = 1; n_minus_half != a_minus_half; ++n_minus_half) {
    term *= xx * kSmallHalfRecips[n_minus_half];
    sum += term;
  }
  // tau_ln + xx guaranteed to be small
  double ee = tt * exp(tau_ln_plus_x) + sum;
  return log(ee) - xx;
}

double ln_sum(double aa, double bb) {
  if (aa > bb) {
    const double tmp = aa;
    aa = bb;
    bb = tmp;
  }
  const double diff = bb - aa;
  if ((aa == -DBL_MAX) || (diff >= kLn2 * 53)) {
    return bb;
  }
  return bb + log(1 + exp(-diff));
}


double beta_small_b_large_a_series_ln(double aa, double bb, double xx, double yy, double s0_ln, uint32_t inv) {
  // BGRAT in DiDonato and Morris.

  // normalized always true, mult always 1
  // a >= 15, b == 0.5 or 1, and ibeta_imp was patched to ensure x > 0.7.
  double bm1 = bb - 1.0;
  double tt = aa + bm1 * 0.5;
  double lx;
  if (yy < 0.35) {
    lx = log1p(-yy);
  } else {
    lx = log(xx);
  }
  double uu = -tt * lx;
  double hh_ln = regularized_gamma_prefix_ln(bb, uu);
  double prefix_ln = hh_ln - log(tgamma_delta_ratio(aa, bb)) - bb * log(tt);
  // validated up to this point

  double pp[15]; // ~8-15 digit accuracy
  pp[0] = 1.0;

  double jj;
  {
    double dxx;
    if (bb == 0.5) {
      // bugfix (17 Jun 2019): original expression could underflow
      // jj = finite_half_gamma_q2(0, uu, nullptr);
      dxx = finite_half_gamma_q2_ln(0, uu);
    } else {
      assert(bb == 1.0);
      dxx = -uu;
    }
    // Underflow is harmless here, jj becomes nonzero for the other terms in
    // the series in that case.
    jj = exp(dxx - hh_ln);
  }
  double sum = jj; // patch in s0 and prefix at the end
  uint32_t tnp1 = 1;
  double lx2 = lx * 0.5;
  lx2 *= lx2;
  double lxp = 1.0;
  double t4_recip = 0.25 / (tt * tt);
  double b2n = bb;
  for (uint32_t nn = 1; nn < 15; ++nn) {
    tnp1 += 2;
    double new_pn = 0.0;
    uint32_t tmp1 = 3;
    const double nn_d = u31tod(nn);
    for (uint32_t mm = 1; mm < nn; ++mm) {
      double mbn = u31tod(mm) * bb - nn_d;
      new_pn += mbn * pp[nn - mm] * kFactorialRecips[tmp1];
      tmp1 += 2;
    }
    new_pn /= nn_d;
    new_pn += bm1 * kFactorialRecips[tnp1];
    pp[nn] = new_pn;

    jj = (b2n * (b2n + 1) * jj + (uu + b2n + 1) * lxp) * t4_recip;
    lxp *= lx2;
    b2n += 2;

    double rr = new_pn * jj;
    sum += rr;
    if (rr > 1) {
      if (fabs(rr) < fabs(kEpsilon * sum)) {
        break;
      }
    } else {
      if (fabs(rr * (1.0 / kEpsilon)) < fabs(sum)) {
        break;
      }
    }
  }
  double result_ln = ln_sum(s0_ln, prefix_ln + log(sum));
  if (!inv) {
    return result_ln;
  }
  return log1p(-exp(result_ln));
}
double ibeta_power_terms_ln(double aa, double bb, double xx, double yy) {
  // returns log((x^a)(y^b) / Beta(a,b))
  //
  // normalized always true
  // prefix always 1
  double cc = aa + bb;

  const double agh = aa + kLanczosG - 0.5;
  const double bgh = bb + kLanczosG - 0.5;
  const double cgh = cc + kLanczosG - 0.5;
  double numer_a;
  double denom_a = lanczos_sum_expg_scaled_imp(aa, &numer_a);
  double numer_b;
  double denom_b = lanczos_sum_expg_scaled_imp(bb, &numer_b);
  double denom_c;
  double numer_c = lanczos_sum_expg_scaled_imp(cc, &denom_c);
  double result = (numer_a * numer_b * numer_c) / (denom_a * denom_b * denom_c);
  result *= sqrt(agh * bgh * kRecipE / cgh);
  double result_ln = log(result);
  double l1 = (xx * bb - yy * agh) / agh;
  double l2 = (yy * aa - xx * bgh) / bgh;
  // Since we're returning log(result) rather than the original result (thus,
  // no intermediate overflow/underflow problems), and we only need 6 digits of
  // precision, we shouldn't need any of the numerous cases in the Boost code.
  return result_ln + aa * log1p(l1) + bb * log1p(l2);
}

double ibeta_a_step_ln(double aa, double bb, double xx, double yy, uint32_t kk) {
  double prefix_ln = ibeta_power_terms_ln(aa, bb, xx, yy) - log(aa);
  double sum = 1.0;
  double term = 1.0;
  const double a_plus_b = aa + bb;
  const double k_minus_1pt5 = u31tod(kk) - 1.5;
  for (double i_d = 0.0; i_d < k_minus_1pt5; ) {
    term *= (a_plus_b + i_d) * xx;
    i_d += 1.0;
    term /= aa + i_d;
    sum += term;
  }
  return prefix_ln + log(sum);
}

double ibeta_fraction2_ln(double aa, double bb, double xx, double yy, uint32_t inv) {
  // normalized always true

  // todo: original DiDonato and Morris paper notes that "x must also be a
  // sufficient distance from p when a > 100"; check if we have a problem there

  double result_ln = ibeta_power_terms_ln(aa, bb, xx, yy);

  // see Boost continued_fraction_b()
  const double ay_minus_bx_plus1 = aa * yy - bb * xx + 1.0;
  double cc = (aa * ay_minus_bx_plus1) / (aa + 1.0);
  if (fabs(cc) < kLentzFpmin) {
    cc = kLentzFpmin;
  }
  const double a_plus_b = aa + bb;
  const double x2 = xx * xx;
  result_ln -= log(cc);
  double dd = 0.0;
  double mm = 0.0;
  while (1) {
    double cur_a = (aa + mm) * (a_plus_b + mm);
    mm += 1.0;
    cur_a *= mm * (bb - mm) * x2;
    double denom = aa + 2 * mm - 1.0;
    cur_a /= denom * denom;
    double cur_b = mm;
    cur_b += (mm * (bb - mm) * xx) / (aa + 2 * mm - 1.0);
    cur_b += ((aa + mm) * (ay_minus_bx_plus1 + mm * (2.0 - xx))) / (aa + 2 * mm + 1.0);
    dd = cur_b + cur_a * dd;
    if (fabs(dd) < kLentzFpmin) {
      dd = kLentzFpmin;
    }
    cc = cur_b + cur_a / cc;
    if (fabs(cc) < kLentzFpmin) {
      cc = kLentzFpmin;
    }
    dd = 1.0 / dd;
    const double delta = cc * dd;
    result_ln -= log(delta);
    if (fabs(delta - 1.0) <= kContinuedFractionEpsilon) {
      return inv? log1p(-exp(result_ln)) : result_ln;
    }
  }
}

double binomial_ccdf_ln(uint32_t nn, uint32_t kk, double xx, double yy, uint32_t inv) {
  // x^n + (n choose 1)x^{n-1}y + ... + (n choose (n-k-1))x^{k+1}y^{n-k-1}
  // This is currently just designed to work with ibeta_imp2_ln().  So we
  // assume x and y positive, y = 1 - x, (n-k) < 40, (n-k) <= (k+1)(y/x).

  // Thanks to the (n-k) <= (k+1)(y/x) condition, it's always reasonable to
  // start from the rightmost term normalized to 1, and sum leftward; then we
  // multiply by the rightmost term at the end.
  const double n_d = u31tod(nn);
  const double k_plus1_d = u31tod(kk + 1);
  const double x_div_y = xx / yy;
  double multiplier = 1.0;
  double cur_term = 1.0;
  // need to shift to avoid overflow for n > 1 billion, (n-k) close to 40
  double shifted_inv_binom_coeff = 1.0 / kExactTestBias;
  double i_d = n_d - k_plus1_d;
  for (uint32_t uii = nn - kk - 1; uii; --uii) {
    //   (n choose (i - 1)) / (n choose i)
    //   i!(n-i)! / (i-1)!(n-i+1)!
    // = i / (n-i+1)
    if (cur_term < kRecip2m53) {
      // Only need to finish computing shifted_inv_binom_coeff: multiply by
      // i!, divide by n!/(n-i)!.
      // Since we support n up to 2^31, n!/(n-i)! may overflow for i>33.
      for (; uii > 33; --uii) {
        shifted_inv_binom_coeff *= i_d;
        i_d -= 1.0;
        shifted_inv_binom_coeff /= n_d - i_d;
      }
      double denom = kFactorialRecips[uii];
      double n_minus_i = n_d - i_d;
      for (; uii; --uii) {
        n_minus_i += 1.0;
        denom *= n_minus_i;
      }
      shifted_inv_binom_coeff /= denom;
      break;
    }
    double cur_binom_ratio = i_d;
    i_d -= 1.0;
    cur_binom_ratio /= n_d - i_d;
    shifted_inv_binom_coeff *= cur_binom_ratio;
    cur_term *= cur_binom_ratio * x_div_y;
    multiplier += cur_term;
  }
  // 1.0 / kExactTestBias = 2^83
  // tried taking log(multiplier / shifted_inv_binom_coeff), but that was
  // slower on my Mac?
  const double result_ln = log(multiplier) + k_plus1_d * log(xx) + (n_d - k_plus1_d) * log(yy) - log(shifted_inv_binom_coeff) + 83 * kLn2;
  if (!inv) {
    return result_ln;
  }
  // The nightmare scenario of tiny result_ln shouldn't happen thanks to the
  // (n-k) <= (k+1)(y/x) condition.
  return log1p(-exp(result_ln));
}


// the '2' in imp2 refers to df1 and df2 arguments instead of a and b
double ibeta_imp2_ln(uint32_t df1, uint32_t df2, double xx, uint32_t inv) {
  // In addition to Boost beta.hpp and its dependencies, see DiDonato and
  // Morris's original paper at
  // https://apps.dtic.mil/dtic/tr/fulltext/u2/a210118.pdf .

  // normalized always true
  //
  // assume a >= 0, b >= 0, a+b > 0, x in [0, 1], a and b are multiples of 0.5
  // for now
  // in practice, x always <= 0.5 right now
  if (df1 == 0) {
    return inv? -DBL_MAX : 0.0;
  }
  if (df2 == 0) {
    return inv? 0.0 : -DBL_MAX;
  }
  if (xx == 0.0) {
    return inv? 0.0 : -DBL_MAX;
  }
  if (xx == 1.0) {
    return inv? -DBL_MAX : 0.0;
  }
  double yy = 1.0 - xx;
  if ((df1 == 1) && (df2 == 1)) {
    return log(asin(sqrt(inv? yy : xx)) * (2.0 / kPi));
  }
  if (df1 == 2) {
    df1 = df2;
    df2 = 2;

    const double tmp = xx;
    xx = yy;
    yy = tmp;

    inv = !inv;
  }
  double aa = u31tod(df1) * 0.5;
  if (df2 == 2) {
    if (df1 == 2) {
      return log(inv? yy : xx);
    }
    double ln_pp;
    if (yy < 0.5) {
      //       pp = inv? -expm1(aa * log1p(-yy)) : exp(aa * log1p(-yy))
      // -> ln_pp = inv? log(-exp(aa * log1p(-yy)) + 1.0) : (aa * log1p(-yy))
      //          = inv? log1p(-exp(aa * log1p(-yy))) : (aa * log1p(-yy))
      const double inner_term = aa * log1p(-yy);
      ln_pp = inv? log1p(-exp(inner_term)) : inner_term;
    } else {
      //       pp = inv? -powm1(xx, aa) : pow(xx, aa);
      // -> ln_pp = inv? log(-powm1(xx, aa)) : (aa * log(xx))
      ln_pp = inv? neg_powm1_imp_ln(xx, aa) : (aa * log(xx));
    }
    return ln_pp;
  }

  double bb = u31tod(df2) * 0.5;
  if ((df1 == 1) || (df2 == 1)) {
    if (xx > 0.5) {
      const uint32_t tmp1 = df1;
      df1 = df2;
      df2 = tmp1;

      double tmp2 = aa;
      aa = bb;
      bb = tmp2;

      tmp2 = xx;
      xx = yy;
      yy = tmp2;

      inv = !inv;
    }
    // can ignore max(a, b) <= 1 branch
    if ((df2 == 1) || ((xx < 0.1) && (bb * xx <= 0.49))) {
      return ibeta_series_ln(aa, bb, xx, inv);
    }
    // a/b swapped, x/y swapped, inv swapped
    if (xx >= 0.3) {
      return ibeta_series_ln(bb, aa, yy, !inv);
    }
    if (bb >= 15) {
      return beta_small_b_large_a_series_ln(bb, aa, yy, xx, -DBL_MAX, !inv);
    }
    double fract_ln = ibeta_a_step_ln(bb, aa, yy, xx, 20);
    return beta_small_b_large_a_series_ln(bb + 20.0, aa, yy, xx, fract_ln, !inv);
  }
  double lambda;
  if (aa < bb) {
    lambda = aa - (aa + bb) * xx;
  } else {
    lambda = (aa + bb) * yy - bb;
  }
  if (lambda < 0.0) {
    const uint32_t tmp1 = df1;
    df1 = df2;
    df2 = tmp1;

    double tmp2 = aa;
    aa = bb;
    bb = tmp2;

    tmp2 = xx;
    xx = yy;
    yy = tmp2;

    inv = !inv;
  }

  // a > 1, b > 1 guaranteed if we get here.
  if (df2 >= 80) {
    return ibeta_fraction2_ln(aa, bb, xx, yy, inv);
  }
  if ((!(df1 % 2)) && (!(df2 % 2)) && (yy != 1.0)) {
    const uint32_t a_int = df1 / 2;
    const uint32_t b_int = df2 / 2;
    const uint32_t kk = a_int - 1;
    const uint32_t nn = b_int + kk;
    return binomial_ccdf_ln(nn, kk, xx, yy, inv);
  }
  // Changed from b*x <= 0.7 since BGRAT has problems with small x, while BPSER
  // can handle this larger b*x case since b<40 is guaranteed.
  if (xx <= 0.7) {
    return ibeta_series_ln(aa, bb, xx, inv);
  }
  const uint32_t nn = (df2 - 1) / 2;
  const double bbar = bb - u31tod(nn);
  double fract_ln = ibeta_a_step_ln(bbar, aa, yy, xx, nn);
  if (df1 > 30) {
    return beta_small_b_large_a_series_ln(aa, bbar, xx, yy, fract_ln, inv);
  }
  fract_ln = ln_sum(fract_ln, ibeta_a_step_ln(aa, bbar, xx, yy, 20));
  return beta_small_b_large_a_series_ln(aa + 20.0, bbar, xx, yy, fract_ln, inv);
}

// ***** begin TstatToLnP *****


double TstatToLnP(double tt, uint32_t df) {
  const double df_d = u31tod(df);
  const double t2 = tt * tt;
  if (df_d > 2 * t2) {
    const double zz = t2 / (df_d + t2);
    return ibeta_imp2_ln(1, df, zz, 1);
  }
  const double zz = df_d / (df_d + t2);
  return ibeta_imp2_ln(df, 1, zz, 0);
}
// ***** end TstatToLnP *****


// ***** begin ChisqToLnP *****

double finite_gamma_q_ln(uint32_t aa, double xx) {
  // a is a positive integer < 30; max(0.6, a-1) < x < 1e10
  // (e^{-x})(1 + x + x^2/2 + x^3/3! + x^4/4! + ... + x^{a-1}/(a-1)!)
  //
  // logarithm:
  // log(1 + x + ... + x^{a-1}/(a-1)!) - x
  // no overflow or underflow danger for main term thanks to bounds
  double sum = 1.0;
  double term = 1.0;
  for (uint32_t nn = 1; nn < aa; ++nn) {//Here,I change nn != aa to nn < aa due to a warning of array out of bounds
    // division is slow enough that the lookup table speeds up this function by
    // >3x
    term *= xx * kSmallRecips[nn];
    sum += term;
  }
  return log(sum) - xx;
}

double lower_gamma_series(double aa, double zz, double init_value) {
  // evaluate init_value + 1 + (z/(a+1)) + (z^2 / ((a+1)(a+2))) +
  //   (z^3 / ((a+1)(a+2)(a+3))) + ...
  // z shouldn't be larger than (a+170), otherwise overflow is possible
  double result = 1;
  double total = init_value + 1;
  do {
    aa += 1.0;
    result *= zz / aa;
    total += result;
  } while (fabs(result) > (kBigEpsilon * kBigEpsilon));
  return total;
}

double upper_gamma_fraction(double a1, double z1) {
  // No overflow issues, this tends to be close to 1.

  // evaluate a_1 / (b_1 + (a_2 / (b_2 + (a_3 / (b_3 + ...)))))
  // see Boost continued_fraction_a(), upper_incomplete_gamma_fract
  double cur_b = z1 - a1 + 3;

  double hh = cur_b;
  const double a0 = a1 - 1.0;
  if (fabs(hh) < kLentzFpmin) {
    hh = kLentzFpmin;
  }
  double cc = hh;
  double dd = 0.0;
  for (double kk = 2.0; kk <= 100.0; kk += 1.0) {
    const double cur_a = kk * (a1 - kk);
    cur_b += 2.0;
    dd = cur_b + cur_a * dd;
    if (fabs(dd) < kLentzFpmin) {
      dd = kLentzFpmin;
    }
    cc = cur_b + cur_a / cc;
    if (fabs(cc) < kLentzFpmin) {
      cc = kLentzFpmin;
    }
    dd = 1.0 / dd;
    const double delta = cc * dd;
    hh *= delta;
    if (fabs(delta - 1.0) < kContinuedFractionEpsilon) {
      break;
    }
  }
  // const double cont_frac = a0 / hh;
  // return 1 / (z1 - a1 + 1 + cont_frac);
  return hh / (hh * (z1 - a1 + 1) + a0);
}

// maximal error of 1.2e-7
double erfc_fast(double zz) {
  const double tt = 1.0 / (1.0 + 0.5 * zz);
  const double tau = tt * exp(((((((((0.17087277 * tt - 0.82215223) * tt + 1.48851587) * tt - 1.13520398) * tt + 0.27886807) * tt - 0.18628806) * tt + 0.09678418) * tt + 0.37409196) * tt + 1.00002368) * tt - 1.26551223 - zz * zz);
  return tau;
}

double igamma_temme_large(double aa, double xx) {
  // 24-bit precision is fine
  const double sigma = (xx - aa) / aa;
  // abs(sigma) < 0.4
  const double phi = -log1pmx(sigma);
  const double sqrt_a = sqrt(aa);
  const double sqrt_phi = sqrt(phi);
  const double yy = aa * phi;
  double zz = kSqrt2 * sqrt_phi;
  if (xx < aa) {
    zz = -zz;
  }
  double workspace[3];
  workspace[0] = (((((kTemmeC0[6] * zz + kTemmeC0[5]) * zz + kTemmeC0[4]) * zz + kTemmeC0[3]) * zz + kTemmeC0[2]) * zz + kTemmeC0[1]) * zz + kTemmeC0[0];
  workspace[1] = (((kTemmeC1[4] * zz + kTemmeC1[3]) * zz + kTemmeC1[2]) * zz + kTemmeC1[1]) * zz + kTemmeC1[0];
  workspace[2] = (kTemmeC2[2] * zz + kTemmeC2[1]) * zz + kTemmeC2[0];
  const double a_recip = 1 / aa;
  double result = (workspace[2] * a_recip + workspace[1]) * a_recip + workspace[0];
  result *= exp(-yy) / ((kSqrt2 * kSqrtPi) * sqrt_a);
  if (xx < aa) {
    result = -result;
  }
  result += erfc_fast(sqrt_a * sqrt_phi) * 0.5;
  return result;
}


// does not guarantee return value <= 0 for now; caller must do that.
double gamma_incomplete_imp2_ln(uint32_t df, double xx) {
  // normalized = true, invert = false
  assert(df);
  assert(xx >= 0.0);
  const double aa = u31tod(df) * 0.5;
  const uint32_t is_small_a = (df < 60) && (aa <= xx + 1) && (xx < 1e10);
  uint32_t is_int = 0;
  uint32_t is_half_int = 0;
  if (is_small_a) {
    is_half_int = df % 2;
    is_int = !is_half_int;
  }
  uint32_t eval_method;
  if (is_int && (xx > 0.6)) {
    eval_method = 0;
  } else if (is_half_int && (xx > 0.2)) {
    eval_method = 1;
  } else if (xx < kSmallEpsilon) {
    // avoid computing log(0)
    // don't need more precision here, 6 digits is enough
    return 0.0;
  } else if (xx < 1.1) {
    eval_method = 2;
  } else {
    const double x_minus_a = xx - aa;
    uint32_t use_temme = 0;
    if (aa > 20) {
      // sigma = abs((x - a) / a);
      // igamma_temme_large() assumes abs(sigma) < 0.95
      if (aa > 200) {
        // abs(sigma) < sqrt(20 / a) < 0.316...
        use_temme = (20 * aa > x_minus_a * x_minus_a);
      } else {
        // abs(sigma) < 0.4
        const double sigma_times_a = fabs(x_minus_a);
        use_temme = (sigma_times_a < 0.4 * aa);
      }
    }
    if (use_temme) {
      eval_method = 5;
    } else {
      // x - (1 / (3 * x)) < a
      // x * x - (1/3) < a * x
      // x * x - a * x < 1/3
      // x * (x - a) < 1/3
      if (xx * x_minus_a < (1.0 / 3.0)) {
        eval_method = 2;
      } else {
        eval_method = 4;
      }
    }
  }
  switch(eval_method) {
  case 0:
    return finite_gamma_q_ln(df / 2, xx);
  case 1:
    return finite_half_gamma_q2_ln(df / 2, xx);
  case 2:
    {
      const double result_ln = regularized_gamma_prefix_ln(aa, xx);
      if (result_ln < kLogMinValue + 22) {
        // init_value overflows.  Not a big deal, this just ends up getting
        // inverted to pval=1.
        // (+22 since aa could theoretically be as large as 2^31.  Todo: find
        // the smallest result_ln value that could result in a nonzero value
        // being returned.)
        return 0.0;
      }
      const double init_value = -aa * exp(-result_ln);
      const double multiplier = -lower_gamma_series(aa, xx, init_value) / aa;
      return result_ln + log(multiplier);
    }
  case 4:
    {
      const double result1_ln = regularized_gamma_prefix_ln(aa, xx);
      const double result2_ln = log(upper_gamma_fraction(aa, xx));
      return result1_ln + result2_ln;
    }
  case 5:
  default:  // silence compiler warning
    {
      // aa large, fabs(xx - aa) relatively small, so no overflow/underflow
      // issues.
      double result = igamma_temme_large(aa, xx);
      if (xx < aa) {
        result = 1.0 - result;
      }
      return log(result);
    }
  }
}


double ChisqToLnP(double chisq, uint32_t df) {
  return MINV(gamma_incomplete_imp2_ln(df, chisq * 0.5), 0.0);
}

// ***** end ChisqToLnP *****

//***********PLINK_FUNCTION*********



//***CONVERT_LN_TO_STRING(by myself)***
string tstat(int n,double t){
  double index_base_e=TstatToLnP(t,n);
  int index_base_10= ceil(index_base_e*log10(kE))-1;
  double coefficient=pow(10, index_base_e*log10(kE)-index_base_10);
  string c_string=to_string(coefficient);
//string e="e";
  string i_string=to_string(index_base_10);
  return c_string + "e" + i_string;
}

string Chisq(double t) {
  double index_base_e=ChisqToLnP(t,1);
  int index_base_10= ceil(index_base_e*log10(kE))-1;
  double coefficient=pow(10, index_base_e*log10(kE)-index_base_10);
  string c_string=to_string(coefficient);
  //string e="e";
  string i_string=to_string(index_base_10);
  return c_string + "e" + i_string;
}

}

#endif
